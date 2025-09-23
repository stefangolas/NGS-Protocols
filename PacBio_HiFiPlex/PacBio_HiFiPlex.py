from pyhamilton import (HamiltonInterface, start_timer)
from pyhamilton.pipetting import (pip_transfer, multi_dispense, double_aspirate_supernatant_96, ethanol_wash, transfer_96, pip_mix, mix_plate,
                                  shear_plate_96, transfer_96, pip_pool)
from pyhamilton.transport import (transport_resource, GripDirection, GrippedResource, GripperParams)
from pyhamilton.devices import (hhs_set_simulation, hhs_create_star_device, hhs_create_usb_device, hhs_set_temp_param,
                               hhs_start_temp_ctrl, hhs_stop_temp_ctrl, hhs_start_shaker, hhs_stop_shaker,
                               initialize_cpac, set_temperature_target_cpac, start_temperature_control_cpac, get_temperature_cpac)
from pyhamilton.consumables import (tracked_volume_aspirate, ReagentTrackedReservoir60mL, ReagentTrackedPlate96, 
                                    ReagentTrackedFalconCarrier24, generate_reagent_summary, ReagentTrackedBulkPlate)
from pyhamilton.ngs import Protocol, LoadingVis
from pyhamilton.resources import TrackedTips, StackedResources, Plate96, Tip96,Waste96, LayoutManager,TipSupportTracker, layout_item

import time
import os


class HHS:
    """Helper class for Hamilton Heater Shaker devices."""
    def __init__(self, node, sequence, lmgr: LayoutManager):
        self.node = node
        self._sequence = sequence
        self.lmgr = lmgr
        self.resource = layout_item(lmgr, Plate96, sequence)
    
    def layout_name(self):
        return self._sequence


# Implement aspiration and dispense first guesses for height
# Implement liquid class guesses

class PacBioHiFiPlexProtocol(Protocol):
    """
    PacBio HiFiPlex library preparation protocol class for Hamilton liquid handler.
    
    This class encapsulates the complete workflow for PacBio HiFiPlex library preparation
    including DNA shearing, cleanup, repair & A-tailing, adapter ligation, and pooling.
    """
    
    def __init__(self, num_samples=96, sample_volume=50, simulation=False, device_simulation=False):
        """
        Initialize the protocol with sample parameters.
        
        Args:
            num_samples (int): Number of samples to process (default: 96)
            sample_volume (int): Volume of each sample in µL (default: 50)
            simulation (bool): Whether to run in simulation mode (default: False)
        """

        self.available_steps = [
            ("Initialize System", "initialize"),
            ("DNA Shearing", "shear_dna"),
            ("Post-Shear Cleanup", "post_shear_cleanup"),
            ("Repair & A-Tailing", "repair_and_a_tailing"),
            ("Adapter Ligation", "adapter_ligation"),
            ("Pooling Ligation", "pooling_ligation")
        ]

        super().__init__()
        self.num_samples = num_samples
        self.sample_volume = sample_volume
        self._reagent_volumes = None
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
        self.lmgr = LayoutManager('PacBio_MultiPlexLibraryPrepDeck_v1.2.lay')

        # Deck Resources
        self.MIDI_OnMagnet = layout_item(self.lmgr, Plate96, 'MIDI_OnMagnet')
        self.MIDI_OffMagnet = layout_item(self.lmgr, Plate96, 'MIDI_Pipette')
        self.LiquidWaste = layout_item(self.lmgr, Waste96, 'LiquidWaste_MPH')
        self.MIDIWaste = layout_item(self.lmgr, Plate96, 'MIDI_Waste')
        self.HSPWaste = layout_item(self.lmgr, Plate96, 'HSP_Waste')
        self.HSP_Plate = layout_item(self.lmgr, Plate96, 'HSP_Pipette')
        self.HSP_Plate_2 = layout_item(self.lmgr, Plate96, 'HSP_Pipette2')
        self.HSP_Park = layout_item(self.lmgr, Plate96, 'Stack_03_0003')

        # Reagent containers
        self.EthanolReservoir = layout_item(self.lmgr, ReagentTrackedBulkPlate, 'RGT_Ethanol')
        self.MagBeads_Container = layout_item(self.lmgr, ReagentTrackedReservoir60mL, 'rgt_cont_60ml_BC_A00_0001')
        self.ElutionBuffer_Container = layout_item(self.lmgr, ReagentTrackedReservoir60mL, 'rgt_cont_60ml_BC_A00_0002')
        self.EDTA = layout_item(self.lmgr, ReagentTrackedReservoir60mL, 'rgt_cont_60ml_BC_A00_0003')
        self.CPAC_Reagent_Plate = layout_item(self.lmgr, ReagentTrackedPlate96, 'CPAC_HSP_0001')
        self.PoolingTubes = layout_item(self.lmgr, ReagentTrackedFalconCarrier24, 'SMP_CAR_24_15x75_A00_0001')


        # Reagent mapping and positions
        self.magbead_positions =  self.MagBeads_Container.assign_reagent_map('MagBeads', range(8))
        self.ethanol_plate = self.EthanolReservoir.assign_reagent_map('Ethanol', range(96))
        self.post_shear_elution_buffer_positions = self.ElutionBuffer_Container.assign_reagent_map('ElutionBuffer', range(8))
        self.ER_Mix_positions = self.CPAC_Reagent_Plate.assign_reagent_map('ER_Mix', [0])
        self.RGT_LigMix_positions = self.CPAC_Reagent_Plate.assign_reagent_map('RGT_LigMix', [1])
        self.EDTA_positions = self.EDTA.assign_reagent_map('EDTA', range(8))

        # Stacked resources
        self.HSP_Stack = StackedResources.from_prefix(
            tracker_id="BioRadHardShell_Stack4",
            prefix="BioRadHardShell_Stack4",
            lmgr=self.lmgr,
            resource_type=Plate96,
            count=5)

        # Tracked tips
        self.tracked_tips_50uL = TrackedTips.from_prefix(
            tracker_id="TIP_50uLF_L",
            volume_capacity=50,
            prefix="TIP_50uLF_L",
            count=8,
            tip_type=Tip96, 
            lmgr=self.lmgr)

        self.tracked_tips_300uL = TrackedTips.from_prefix(
            tracker_id="STF_L",
            volume_capacity=300,
            prefix="STF_L",
            count=8,
            tip_type=Tip96,
            lmgr=self.lmgr)

        self.tracked_tips_1000uL = TrackedTips.from_prefix(
            tracker_id="HTF_L",
            volume_capacity=1000,
            prefix="HTF_L",
            count=2,
            tip_type=Tip96,
            lmgr=self.lmgr)
        
        self.tip_support_resource = layout_item(self.lmgr, Tip96, 'TipSupport_0001')
        self.tip_support = TipSupportTracker(self.tip_support_resource)

        # Batch together tracked objects for resource consumption estimates
        self.tracked_reagent_vessels = [self.MagBeads_Container, self.ElutionBuffer_Container, self.CPAC_Reagent_Plate, self.EDTA, self.PoolingTubes, self.EthanolReservoir]
        self.tracked_tips = [self.tracked_tips_50uL, self.tracked_tips_300uL, self.tracked_tips_1000uL]
        self.stacked_resources = [self.HSP_Stack]

    
    def _setup_volumes(self):
        """Setup volume parameters based on sample volume."""
        self.magbead_mix_volume = 1000
        self.post_shear_magbead_volume = self.sample_volume
        self.first_supernatant_removal_volume = self.sample_volume * 2
        self.supernatant_removal_volume = self.sample_volume + self.post_shear_magbead_volume
        self.m1_mix_volume = min(self.sample_volume * 1.6, 1000)
        self.post_shear_etoh_wash_volume = 200
        self.post_shear_elution_buffer_volume = 30
        self.post_shear_elution_volume = 25.5

    
    def _setup_hhs_devices(self):
        """Setup HHS (Hamilton Heater Shaker) devices."""
        self.HHS4_MIDI = HHS(node=4, sequence="HHS4_MIDI", lmgr=self.lmgr)
        self.HHS5_MIDI = HHS(node=5, sequence="HHS5_MIDI", lmgr=self.lmgr)
        self.HHS3_MIDI = HHS(node=3, sequence="HHS3_MIDI", lmgr=self.lmgr)
        self.HHS1_HSP = HHS(node=1, sequence="HHS1_HSP", lmgr=self.lmgr)
        self.HHS2_HSP = HHS(node=2, sequence="HHS2_HSP", lmgr=self.lmgr)

    
    def initialize_hhs(self, ham_int):
        """Initialize Hamilton Heater Shaker devices."""
        hhs_set_simulation(ham_int, 0 if not self.device_simulation else 1)
        
        for node in range(1, 6):
            try:
                hhs_create_usb_device(ham_int, node)
                print(f"Created USB device for ML_STAR {node}")
            except Exception as e:
                print(f"Warning: Could not initialize HHS node {node}: {e}")
    
    def initialize(self):
        """Initialize the Hamilton system and all devices."""
        with HamiltonInterface(simulating=self.simulation, server_mode=False, windowed=self.windowed, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            self.initialize_hhs(ham_int)
            initialize_cpac(ham_int, controller_id=1, simulating=False)
            
            self._initialized = True
            print("Protocol initialization completed successfully.")
    
    def shear_dna(self):
        """Step 1: DNA Shearing."""
        if not self._initialized:
            print("Warning: Protocol not initialized. Call initialize() first.")
        
        print("Starting DNA Shearing step...")

        with HamiltonInterface(simulating=self.simulation, server_mode=False, windowed=self.windowed, persistent=self.persistent) as ham_int:
            ham_int.initialize()

            mix_plate(ham_int, self.tracked_tips_300uL, self.tip_support, num_samples=self.num_samples, plate=self.MIDI_OnMagnet,
                      mixing_volume=100, mix_cycles=5, liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
            
            print("DNA Shearing completed.")
    
    def post_shear_cleanup(self, starting_from=0):
        """Step 2: Post-shear magnetic bead cleanup."""
        print("Starting Post-Shear Cleanup step...")

        with HamiltonInterface(simulating=self.simulation, server_mode=False, windowed=self.windowed, persistent=self.persistent) as ham_int:
            ham_int.initialize()

            # Mix magnetic beads
            print("Mixing magnetic beads...")
            pip_mix(ham_int, self.tracked_tips_1000uL, self.magbead_positions, 
                   liquid_class='PacBio_HighVolume_PurificationBeadsMix_SurfaceEmptyV2', 
                   mix_volume=self.magbead_mix_volume, mix_cycles=20, liquid_height=0)

            # Transfer beads to MIDI Off Magnet plate
            print("Transferring beads to MIDI Off Magnet plate...")
            sample_poss = [(self.MIDI_OffMagnet, i) for i in range(self.num_samples)]
            vols = [self.sample_volume] * self.num_samples
            pip_transfer(ham_int, self.tracked_tips_300uL, self.magbead_positions, sample_poss, vols, dispense_height=1,
                         liquid_class='PacBio_SVF_ProNexPurificationBeads_SurfaceEmpty_v3')

            # Mix plate with 96 channel head
            print("Mixing plate with 96 channel head...")
            mix_plate(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples, self.MIDI_OffMagnet, self.m1_mix_volume, 
                     mix_cycles=20, liquid_class='PacBio_SVF_ProNexPurificationBeads_SurfaceEmpty_v3')

            # Start bead incubation timer
            bead_incubation_timer = start_timer(600)

            # Pre-stamp elution buffer
            hsp_plate_positions = [(self.HSP_Plate, i) for i in range(self.num_samples)]
            post_shear_elution_buffer_volumes = [self.post_shear_elution_volume] * self.num_samples
            multi_dispense(ham_int, self.tracked_tips_300uL, self.post_shear_elution_buffer_positions, 
                          hsp_plate_positions, post_shear_elution_buffer_volumes, aspiration_height=1, dispense_height=1,
                          liquid_class='PacBio_SVF_EB_SingleDispense_JetEmpty')

            # Wait for bead incubation
            bead_incubation_timer.wait(skip=True)
            
            # Move MIDI plate to magnet
            transport_resource(ham_int, self.MIDI_OffMagnet, self.MIDI_OnMagnet, 
                             resource_type=GrippedResource.MIDI, core_gripper=True)

            # Bead settling
            bead_settle_timer = start_timer(300)
            bead_settle_timer.wait(skip=True)

            # Mix magnetic beads
            mix_plate(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples, self.MIDI_OnMagnet, self.m1_mix_volume, 
                     liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty', 
                     mix_cycles=6, liquid_height=0)

            # Additional incubation
            additional_incubation_timer = start_timer(180)
            additional_incubation_timer.wait(skip=True)

            # Remove supernatant
            if self.sample_volume > 130:
                transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples, self.MIDI_OnMagnet, self.LiquidWaste, 
                           'StandardVolumeFilter_Water_DispenseSurface_Empty', 
                           self.first_supernatant_removal_volume, aspiration_height=0, dispense_height=10)

            double_aspirate_supernatant_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples, self.MIDI_OnMagnet, 
                                          self.LiquidWaste, first_volume=270, second_volume=30, 
                                          second_aspiration_height=0, liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            # Ethanol wash
            ethanol_wash(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                        ethanol_plate=self.EthanolReservoir, magnet_plate=self.MIDI_OnMagnet, waste_plate=self.LiquidWaste,
                        wash_volume=200, first_removal_volume=242, second_removal_volume=58,
                        liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            # Reset tips
            self.tracked_tips_300uL.reset_all()

            # Air dry ethanol
            air_dry_timer = start_timer(60)
            air_dry_timer.wait(skip=True)

            # Move MIDI plate off magnet
            transport_resource(ham_int, self.MIDI_OnMagnet, self.MIDI_OffMagnet, 
                             resource_type=GrippedResource.MIDI, core_gripper=True)

            # Move pre-stamped elution buffer from HSP plate to MIDI plate
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples, 
                        self.HSP_Plate, self.MIDI_OffMagnet, liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty', volume=26.0, 
                       aspiration_height=0.2, dispense_height=2.0)

            # Move MIDI plate to HHS and shake
            transport_resource(ham_int, self.MIDI_OffMagnet, self.HHS3_MIDI.resource, 
                             resource_type=GrippedResource.MIDI, core_gripper=True)

            # Incubate
            incubate_timer = start_timer(300)
            incubate_timer.wait(skip=True)

            # Plate management
            transport_resource(ham_int, self.HSP_Plate, self.HSPWaste, 
                             resource_type=GrippedResource.PCR, core_gripper=True)
            transport_resource(ham_int, self.HSP_Stack.fetch_next(), self.HSP_Plate, 
                             resource_type=GrippedResource.PCR, core_gripper=True)

            # Additional incubation
            if not self.device_simulation:
                time.sleep(30)

            # Final transfers
            transport_resource(ham_int, self.HHS3_MIDI.resource, self.MIDI_OnMagnet, 
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            
            if not self.device_simulation:
                time.sleep(200)

            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples, 
                        self.MIDI_OnMagnet, self.HSP_Plate, liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                       aspiration_height=0.3, dispense_height=2.0, volume=30)

            transport_resource(ham_int, self.MIDI_OnMagnet, self.MIDIWaste, 
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            

            print("Post-Shear Cleanup completed.")
    
    def repair_and_a_tailing(self, starting_from=0):
        """Step 3: DNA Repair and A-Tailing."""
        print("Starting Repair and A-Tailing step...")

        with HamiltonInterface(simulating=self.simulation, server_mode=False, windowed=self.windowed, persistent=self.persistent) as ham_int:
            ham_int.initialize()

            # Check current temperature
            temperature = get_temperature_cpac(ham_int, controller_id=1, device_id=1)

            # Transfer ER Mix to HSP plate
            HSP_Plate_positions = [(self.HSP_Plate, i) for i in range(self.num_samples)]
            pip_transfer(ham_int, self.tracked_tips_50uL, self.ER_Mix_positions, HSP_Plate_positions, 
                        volumes=[5.5] * self.num_samples, 
                        liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')

            # Mix HSP plate
            mix_plate(ham_int, self.tracked_tips_50uL, self.tip_support, self.num_samples,
                    plate=self.HSP_Plate, mixing_volume=24, liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty',
                    liquid_height=0.0, mix_cycles=30)

            # Set temperature controls
            hhs_start_temp_ctrl(ham_int, device_number=self.HHS1_HSP.node, temperature=37, wait_for_temp_reached=False)
            hhs_start_temp_ctrl(ham_int, device_number=self.HHS2_HSP.node, temperature=65, wait_for_temp_reached=False)

            # Move HSP plate to HHS1 (37°C)
            transport_resource(ham_int, self.HSP_Plate, self.HHS1_HSP.resource, 
                             resource_type=GrippedResource.MIDI, core_gripper=True)

            # First incubation at 37°C
            er_mix_incubate_timer = start_timer(1800)
            er_mix_incubate_timer.wait(skip=self.device_simulation)

            # Stop temperature control on HHS1
            hhs_stop_temp_ctrl(ham_int, device_number=self.HHS1_HSP.node)

            # Move plate from HHS1 to HHS2 (65°C)
            transport_resource(ham_int, self.HHS1_HSP.resource, self.HHS2_HSP.resource, 
                             resource_type=GrippedResource.MIDI, core_gripper=True)

            # Second incubation at 65°C
            er_second_incubate_timer = start_timer(300)
            er_second_incubate_timer.wait(skip=self.device_simulation)

            # Stop temperature control on HHS2
            hhs_stop_temp_ctrl(ham_int, device_number=self.HHS2_HSP.node)

            # Move plate back to HSP_Plate position
            transport_resource(ham_int, self.HHS2_HSP.resource, self.HSP_Plate, 
                             resource_type=GrippedResource.PCR, core_gripper=True)

            print("Repair and A-Tailing completed.")
    
    def adapter_ligation(self, starting_from=0):
        """Step 4: Adapter Ligation."""
        print("Starting Adapter Ligation step...")

        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, persistent=self.persistent) as ham_int:
            ham_int.initialize()

            # Transfer M2 adapters to HSP plate
            HSP_Plate_positions = [(self.HSP_Plate, i) for i in range(self.num_samples)]

            # Transport from HSP stack to HSP_Plate_2
            transport_resource(ham_int, self.HSP_Plate_2, self.HSP_Park, 
                             resource_type=GrippedResource.PCR, core_gripper=True)

            # Transfer Ligation Mix to HSP_Plate_2
            pip_transfer(ham_int, self.tracked_tips_50uL, self.RGT_LigMix_positions, HSP_Plate_positions, 
                        volumes=[10.5] * self.num_samples,
                        liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')

            # Mix HSP plate
            mix_plate(ham_int, self.tracked_tips_50uL, self.tip_support, self.num_samples,
                    self.HSP_Plate, liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty', 
                     mix_cycles=20, mixing_volume=34, liquid_height=0.1)

            # Move HSP plate to HSP_Plate2
            transport_resource(ham_int, self.HSP_Plate, self.HSP_Plate_2, 
                             resource_type=GrippedResource.PCR, core_gripper=True)

            # Ligation incubation (could be done in parallel with other operations)
            if not self.device_simulation:
                time.sleep(1800)  # 30 minutes


            # Multidispense EDTA to HSP_Plate_2 to prepare for MPH stamp
            HSP_Plate_2_positions = [(self.HSP_Plate_2, i) for i in range(self.num_samples)]
            multi_dispense(ham_int, self.tracked_tips_50uL, self.EDTA_positions, HSP_Plate_2_positions, volumes=[5.0] * self.num_samples,
                           liquid_class='Tip_50ulFilter_Water_DispenseJet_Empty')

            # Transfer from HSP_Plate2 to HSP_Plate
            transfer_96(ham_int, self.tracked_tips_50uL, self.tip_support, self.num_samples, self.HSP_Plate_2, self.HSP_Plate, 
                       volume = 50, liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')

            # Mix HSP plate
            mix_plate(ham_int, self.tracked_tips_50uL, self.tip_support, self.num_samples, self.HSP_Plate, 
                     mixing_volume=34, liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty', mix_cycles=3)

            # Move HSP_Plate2 to waste
            transport_resource(ham_int, self.HSP_Plate_2, self.HSPWaste, 
                             resource_type=GrippedResource.PCR, core_gripper=True)

            print("Adapter Ligation completed.")
    
    def pooling_ligation(self, starting_from=0):
        """Step 5: Pool ligated samples."""
        print("Starting Pooling Ligation step...")
        
        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed) as ham_int:
            ham_int.initialize()
            
            HSP_Plate_positions = [(self.HSP_Plate, i) for i in range(self.num_samples)]
            pooling_tube = [(self.PoolingTubes, 1)]
            
            # Transfer from HSP_Plate to PoolingTubes
            pip_pool(ham_int, self.tracked_tips_50uL, HSP_Plate_positions, pooling_tube, volumes=[50] * self.num_samples,
                       liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')

            print("Pooling Ligation completed.")
    
    def run_complete_protocol(self):
        """Run the complete PacBio HiFiPlex protocol."""
        print(f"Starting PacBio HiFiPlex protocol for {self.num_samples} samples...")
        print(f"Sample volume: {self.sample_volume} µL")
        print(f"Simulation mode: {self.simulation}")
        
        try:
            # Run protocol steps
            self.initialize()
            self.shear_dna()
            self.post_shear_cleanup()
            self.repair_and_a_tailing()
            self.adapter_ligation()
            self.pooling_ligation()
            
            print("PacBio HiFiPlex protocol completed successfully!")
            
        except Exception as e:
            print(f"Protocol failed with error: {e}")
            raise
    

# Example usage
if __name__ == "__main__":
    
    # Create protocol instance
    protocol = PacBioHiFiPlexProtocol(device_simulation=True)
    protocol.run_protocol()

