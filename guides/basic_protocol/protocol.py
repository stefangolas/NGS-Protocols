"""Simple NGS Protocol - Protocol Class Implementation"""

import os
from pathlib import Path
from pyhamilton import (HamiltonInterface, LayoutManager, layout_item, 
                        start_timer, normal_logging)
from pyhamilton.pipetting import (pip_transfer, transfer_96, pip_mix, 
                                  double_aspirate_supernatant_96, ethanol_wash)
from pyhamilton.consumables import (ReagentTrackedReservoir60mL, 
                                    ReagentTrackedBulkPlate, 
                                    ReagentTrackedEppiCarrier32)
from pyhamilton.resources import (Plate96, Tip96, TrackedTips, 
                                  TipSupportTracker, StackedResources, Lid)
from pyhamilton.transport import transport_resource, GrippedResource, GripDirection
from pyhamilton.devices import (hhs_set_simulation, hhs_create_usb_device, 
                                hhs_start_shaker, hhs_stop_shaker,
                                odtc_connect, odtc_initialize, odtc_open_door,
                                odtc_close_door, odtc_execute_protocol, 
                                odtc_wait_for_idle)
from pyhamilton.ngs import Protocol, generate_tadm_report


def get_parent_lay_file():
    """Find the .lay file in the parent directory."""
    parent = Path(__file__).resolve().parent.parent
    for file in parent.glob("*.lay"):
        return str(file)
    return None


class HHS:
    """Hamilton Heater Shaker wrapper class."""
    def __init__(self, node, sequence, lmgr):
        self.node = node
        self._sequence = sequence
        self.lmgr = lmgr
        self.resource = layout_item(lmgr, Plate96, sequence)
    
    def layout_name(self):
        return self._sequence


class BasicProtocol(Protocol):
    """
    Simple NGS protocol demonstrating basic workflow:
    - Reagent dispensing
    - Thermal cycling
    - Magnetic bead cleanup with ethanol washes
    """
    
    def __init__(self, num_samples=8, sample_volume=50, simulation=False, device_simulation=False):
        """
        Initialize the protocol with sample parameters.
        
        Args:
            num_samples (int): Number of samples to process (default: 8)
            sample_volume (int): Volume of each sample in µL (default: 50)
            simulation (bool): Whether to run in simulation mode (default: False)
            device_simulation (bool): Whether to simulate devices (default: False)
        """
        
        self.available_steps = [
            ("Initialize System", "initialize"),
            ("Add Enzyme Mix", "add_enzyme_mix"),
            ("Thermal Cycling", "thermal_cycling"),
            ("Magnetic Bead Cleanup", "magnetic_bead_cleanup")
        ]
        
        super().__init__()
        self.num_samples = num_samples
        self.sample_volume = sample_volume
        self.simulation = simulation
        self.device_simulation = device_simulation
        
        # Initialize deck layout and resources
        self._setup_deck_layout()
        self._setup_volumes()
        self._setup_hhs_devices()
        self._initialized = False
    
    def _setup_deck_layout(self):
        """Setup deck layout and resources."""
        # Deck layout import
        lay_file = get_parent_lay_file()
        self.lmgr = LayoutManager(lay_file)
        
        # Deck Resources - Plates
        self.HSP_CPAC = layout_item(self.lmgr, Plate96, 'HSP_CPAC')
        self.HSP_Pipette2 = layout_item(self.lmgr, Plate96, 'HSP_Pipette2')
        self.HSP_ODTC = layout_item(self.lmgr, Plate96, 'HSP_ODTC')
        self.HSP_ODTC_Lid = layout_item(self.lmgr, Lid, 'Ham_ComfortLid_ODTC')
        
        # MIDI Plates
        self.MIDI_OnMagnet = layout_item(self.lmgr, Plate96, 'MIDI_OnMagnet')
        self.MIDI_Waste = layout_item(self.lmgr, Plate96, 'MIDI_Waste')
        
        # Sample Plate
        self.HSP_Pipette2.assign_label('Samples')
        
        # Small volume reagents - Use EppiCarrier32
        self.CAR_VIALS_SMALL = layout_item(self.lmgr, ReagentTrackedEppiCarrier32, 'CAR_VIALS_SMALL')
        self.Enzyme_Mix_position = self.CAR_VIALS_SMALL.assign_reagent_map('Enzyme_Mix', [1])
        
        # Bulk reagents - Use reservoirs
        self.RGT_01 = layout_item(self.lmgr, ReagentTrackedReservoir60mL, 'RGT_01')
        self.RGT_02 = layout_item(self.lmgr, ReagentTrackedReservoir60mL, 'RGT_02')
        self.Magnetic_Beads_positions = self.RGT_01.assign_reagent_map('Magnetic_Beads', range(8))
        self.Elution_Buffer_positions = self.RGT_02.assign_reagent_map('Elution_Buffer', range(8))
        
        self.EthanolReservoir = layout_item(self.lmgr, ReagentTrackedBulkPlate, 'Ethanol_Reservoir')
        
        # Stacked resources
        self.MIDI_Stack = StackedResources.from_prefix(
            tracker_id="ABGENE_MIDI_Stack1", 
            prefix="ABGENE_MIDI_Stack1",
            count=3, 
            lmgr=self.lmgr, 
            resource_type=Plate96
        )
        
        self.Lid_Stack = StackedResources.from_prefix(
            tracker_id="Ham_ComfortLid_Stack_ParkPos", 
            prefix="Ham_ComfortLid_Stack_ParkPos",
            count=4, 
            lmgr=self.lmgr, 
            resource_type=Lid
        )
        
        # Tracked tips
        self.tracked_tips_50uL = TrackedTips.from_prefix(
            tracker_id="TIP_50ulF_L",
            volume_capacity=50,
            prefix="TIP_50ulF_L",
            count=8,
            tip_type=Tip96,
            lmgr=self.lmgr
        )
        
        self.tracked_tips_300uL = TrackedTips.from_prefix(
            tracker_id="stf_l",
            volume_capacity=300,
            prefix="stf_l",
            count=8,
            tip_type=Tip96,
            lmgr=self.lmgr
        )
        
        # Tip support tracker
        self.tip_support_resource = layout_item(self.lmgr, Tip96, 'TipSupport_0001')
        self.tip_support = TipSupportTracker(self.tip_support_resource)
        
        # Batch tracked objects for resource consumption
        self.tracked_reagent_vessels = [self.CAR_VIALS_SMALL, self.RGT_01, self.HSP_Pipette2]
        self.tracked_tips = [self.tracked_tips_50uL, self.tracked_tips_300uL]
        self.stacked_resources = [self.MIDI_Stack, self.Lid_Stack]
    
    def _setup_volumes(self):
        """Setup volume parameters based on sample volume."""
        self.enzyme_mix_volume = 10
        self.bead_volume = 60
        self.ethanol_wash_volume = 200
        self.elution_volume = 30
    
    def _setup_hhs_devices(self):
        """Setup HHS (Hamilton Heater Shaker) devices."""
        self.HHS3_MIDI = HHS(node=3, sequence="HHS3_MIDI", lmgr=self.lmgr)
    
    def initialize(self):
        """Initialize the Hamilton system and all devices."""
        print("Starting system initialization...")
        
        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, 
                             server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())
            
            # Initialize HHS devices
            print("Initializing HHS devices...")
            hhs_set_simulation(ham_int, 1 if self.device_simulation else 0)
            try:
                hhs_create_usb_device(ham_int, 3)
                print("  HHS node 3: OK")
            except Exception as e:
                print(f"  Warning: Could not initialize HHS node 3: {e}")
            
            # Initialize ODTC
            print("Initializing thermal cycler...")
            self.device_id = odtc_connect(ham_int, simulation_mode=self.device_simulation,
                                         local_ip='1.2.3.4', device_ip='5.6.7.8', 
                                         device_port='COM4')
            odtc_initialize(ham_int, device_id=self.device_id)
            
            print(f"Processing {self.num_samples} samples")
            
            self._initialized = True
            print("Protocol initialization completed successfully.\n")
    
    def add_enzyme_mix(self):
        """Step 1: Add enzyme mix to samples."""
        print("=== Step 1: Adding Enzyme Mix ===\n")
        
        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, 
                             server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())
            
            # Move sample plate to CPAC for reagent addition
            print("Moving sample plate to CPAC position...")
            transport_resource(ham_int, self.HSP_Pipette2, self.HSP_CPAC,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Prepare sample positions
            HSP_CPAC_positions = [(self.HSP_CPAC, idx) for idx in range(self.num_samples)]
            
            # Add enzyme mix to each sample
            print(f"Dispensing {self.enzyme_mix_volume} µL enzyme mix to {self.num_samples} samples...")
            pip_transfer(ham_int, self.tracked_tips_50uL, self.Enzyme_Mix_position,
                        HSP_CPAC_positions, [self.enzyme_mix_volume] * self.num_samples,
                        mix_cycles=5, vol_mix_dispense=10,
                        liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty',
                        aspiration_height=1, dispense_height=1)
            
            # Move plate back to pipetting position
            print("Returning plate to pipetting position...")
            transport_resource(ham_int, self.HSP_CPAC, self.HSP_Pipette2,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            print("Enzyme mix addition completed.\n")
    
    def thermal_cycling(self):
        """Step 2: Run thermal cycling protocol."""
        print("=== Step 2: Thermal Cycling ===\n")
        
        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, 
                             server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())
            
            # Load plate to thermal cycler
            print("Loading plate to thermal cycler...")
            odtc_open_door(ham_int, device_id=self.device_id)
            transport_resource(ham_int, self.HSP_Pipette2, self.HSP_ODTC,
                             grip_direction=GripDirection.RIGHT,
                             resource_type=GrippedResource.PCR, iswap=True)
            
            # Add lid
            print("Adding thermal cycler lid...")
            transport_resource(ham_int, self.Lid_Stack.fetch_next(), self.HSP_ODTC_Lid,
                             grip_direction=GripDirection.RIGHT,
                             resource_type=GrippedResource.LID, iswap=True)
            
            odtc_close_door(ham_int, device_id=self.device_id)
            
            # Run thermal cycling program
            print("Running thermal cycling program...")
            odtc_execute_protocol(ham_int, device_id=self.device_id,
                                method_name='EnzymaticReaction.xml', 
                                simulating=self.device_simulation)
            
            print("Waiting for thermal cycling to complete...")
            odtc_wait_for_idle(ham_int, device_id=self.device_id, 
                              simulating=self.device_simulation, check_interval=5)
            
            # Unload plate from thermal cycler
            print("Unloading plate from thermal cycler...")
            odtc_open_door(ham_int, device_id=self.device_id)
            transport_resource(ham_int, self.HSP_ODTC_Lid, self.Lid_Stack.put_back(),
                             grip_direction=GripDirection.RIGHT,
                             resource_type=GrippedResource.LID, iswap=True)
            transport_resource(ham_int, self.HSP_ODTC, self.HSP_Pipette2,
                             grip_direction=GripDirection.RIGHT,
                             resource_type=GrippedResource.PCR, iswap=True)
            odtc_close_door(ham_int, device_id=self.device_id)
            
            print("Thermal cycling completed.\n")
    
    def magnetic_bead_cleanup(self):
        """Step 3: Perform magnetic bead cleanup with ethanol washes."""
        print("=== Step 3: Magnetic Bead Cleanup ===\n")
        
        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, 
                             server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())
            
            # Get fresh MIDI plate and add beads
            print("Preparing MIDI plate with magnetic beads...")
            transport_resource(ham_int, self.MIDI_Stack.fetch_next(), self.HHS3_MIDI.resource,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            
            HHS3_MIDI_positions = [(self.HHS3_MIDI.resource, idx) for idx in range(self.num_samples)]
            
            # Mix beads before transferring
            print("Mixing magnetic beads...")
            pip_mix(ham_int, tips=self.tracked_tips_300uL, 
                    positions_to_mix=self.Magnetic_Beads_positions,
                    mix_volume=50, mix_cycles=10,
                    liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                    liquid_height=1)
            
            # Add beads to MIDI plate
            print(f"Adding {self.bead_volume} µL beads to MIDI plate...")
            pip_transfer(ham_int, self.tracked_tips_300uL, self.Magnetic_Beads_positions,
                        HHS3_MIDI_positions, [self.bead_volume] * self.num_samples,
                        liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                        aspiration_height=1, dispense_height=1)
            
            # Transfer samples to beads
            print(f"Transferring {self.sample_volume} µL samples to beads...")
            transfer_96(ham_int, self.tracked_tips_300uL, tip_support=self.tip_support,
                       num_samples=self.num_samples, source_plate=self.HSP_Pipette2,
                       target_plate=self.HHS3_MIDI.resource, volume=self.sample_volume,
                       liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                       aspiration_height=1, dispense_height=1)
            
            # Mix samples with beads
            print("Mixing samples with beads (10 sec @ 1000 rpm)...")
            hhs_start_shaker(ham_int, self.HHS3_MIDI.node, 1000)
            shake_timer = start_timer(10)
            shake_timer.wait(skip=self.device_simulation)
            hhs_stop_shaker(ham_int, self.HHS3_MIDI.node)
            
            # Move to magnet
            print("Moving plate to magnetic stand...")
            transport_resource(ham_int, self.HHS3_MIDI.resource, self.MIDI_OnMagnet,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Let beads settle
            print("Settling beads on magnet (60 sec)...")
            settle_timer = start_timer(60)
            settle_timer.wait(skip=self.device_simulation)
            
            # Remove supernatant
            print("Removing supernatant...")
            double_aspirate_supernatant_96(ham_int, self.tracked_tips_300uL, self.tip_support,
                                          self.num_samples, source_plate=self.MIDI_OnMagnet,
                                          destination_plate=self.MIDI_Waste,
                                          first_volume=100, second_volume=20,
                                          liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                                          second_aspiration_height=0.8, dispense_height=5)
            
            # Ethanol washes
            print("Performing ethanol washes (2x)...")
            for wash_num in range(2):
                print(f"  Ethanol wash {wash_num + 1}/2...")
                ethanol_wash(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                            ethanol_plate=self.EthanolReservoir, magnet_plate=self.MIDI_OnMagnet,
                            waste_plate=self.MIDI_Waste, wash_volume=self.ethanol_wash_volume,
                            first_removal_volume=180, second_removal_volume=30,
                            liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')
            
            # Air dry
            print("Air drying beads (60 sec)...")
            dry_timer = start_timer(60)
            dry_timer.wait(skip=self.device_simulation)
            
            # Add elution buffer
            print(f"Adding {self.elution_volume} µL elution buffer...")
            MIDI_OnMagnet_positions = [(self.MIDI_OnMagnet, idx) for idx in range(self.num_samples)]
            pip_transfer(ham_int, self.tracked_tips_300uL, self.Elution_Buffer_positions,
                        MIDI_OnMagnet_positions, [self.elution_volume] * self.num_samples,
                        liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                        aspiration_height=1, dispense_height=1)
            
            # Mix for elution
            print("Mixing for elution (30 sec @ 1000 rpm)...")
            transport_resource(ham_int, self.MIDI_OnMagnet, self.HHS3_MIDI.resource,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            
            hhs_start_shaker(ham_int, self.HHS3_MIDI.node, 1000)
            shake_timer = start_timer(30)
            shake_timer.wait(skip=self.device_simulation)
            hhs_stop_shaker(ham_int, self.HHS3_MIDI.node)
            
            # Return to magnet
            print("Moving to magnet for final collection...")
            transport_resource(ham_int, self.HHS3_MIDI.resource, self.MIDI_OnMagnet,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            
            print("Settling beads (60 sec)...")
            settle_timer = start_timer(60)
            settle_timer.wait(skip=self.device_simulation)
            
            # Transfer purified eluate to output plate
            print("Collecting purified product...")
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                       source_plate=self.MIDI_OnMagnet, target_plate=self.HSP_Pipette2,
                       volume=self.elution_volume - 5,
                       liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                       aspiration_height=0.7, dispense_height=1)
            
            print("Magnetic bead cleanup completed.\n")
    
    def run_complete_protocol(self):
        """Run the complete simple NGS protocol."""
        print(f"\n{'='*60}")
        print(f"Starting Simple NGS Protocol")
        print(f"{'='*60}")
        print(f"Number of samples: {self.num_samples}")
        print(f"Sample volume: {self.sample_volume} µL")
        print(f"Simulation mode: {self.simulation}")
        print(f"Device simulation: {self.device_simulation}")
        print(f"{'='*60}\n")
        
        try:
            # Run protocol steps
            self.initialize()
            self.add_enzyme_mix()
            self.thermal_cycling()
            self.magnetic_bead_cleanup()
            
            print(f"\n{'='*60}")
            print("Simple NGS Protocol completed successfully!")
            print(f"Purified product in HSP_Pipette2")
            print(f"Volume: {self.elution_volume - 5} µL per sample")
            print(f"{'='*60}\n")
            
        except Exception as e:
            print(f"\nProtocol failed with error: {e}")
            raise


# Example usage
if __name__ == "__main__":
    # Create protocol instance
    protocol = BasicProtocol(
        num_samples=8,
        sample_volume=50,
        simulation=True,
        device_simulation=True
    )
    
    # Run the complete protocol
    protocol.run_protocol()
    
    # Generate TADM report
    generate_tadm_report()