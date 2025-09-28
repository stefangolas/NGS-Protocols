from pyhamilton import (HamiltonInterface, LayoutManager, Plate96, Tip96, hhs_set_simulation, normal_logging, start_timer)

from pyhamilton.devices import (initialize_cpac, set_temperature_target_cpac, start_temperature_control_cpac,
                                 get_temperature_cpac, odtc_initialize, odtc_connect, odtc_open_door, odtc_wait_for_idle,
                                 odtc_close_door, odtc_execute_protocol,
                                 hhs_create_usb_device, hhs_start_shaker, hhs_stop_shaker)

from pyhamilton.consumables import (tracked_volume_aspirate, ReagentTrackedReservoir60mL, ReagentTrackedPlate96, 
                                    ReagentTrackedFalconCarrier24, generate_reagent_summary, ReagentTrackedBulkPlate, ReagentTrackedEppiCarrier32)
from pyhamilton.ngs import Protocol, LoadingVis, generate_tadm_report

from pyhamilton.resources import TrackedTips, StackedResources, Plate96, Tip96,Waste96, LayoutManager,TipSupportTracker, layout_item, Lid

from pyhamilton.pipetting import (pip_transfer, multi_dispense, double_aspirate_supernatant_96, ethanol_wash, transfer_96, pip_mix, mix_plate,
                                  shear_plate_96, transfer_96)
from pyhamilton.transport import (transport_resource, GripDirection, GrippedResource, GripperParams)

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


class QIAseqRNAFusionProtocol(Protocol):
    """
    QIAseq RNA Fusion XP Panels library preparation protocol for Hamilton liquid handler.
    
    This class implements the complete workflow for QIAseq RNA Fusion library preparation
    including First Strand DNA Synthesis, Second Strand DNA Synthesis, End Repair & A-Tailing,
    Adapter Ligation, Sample Cleanup, SPE Target Enrichment, and Universal PCR.
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
            ("First Strand DNA Synthesis", "first_strand_dna_synthesis"),
            ("Second Strand DNA Synthesis", "second_strand_dna_synthesis"),
            ("End Repair & A-Tailing", "end_repair_a_tailing"),
            ("Adapter Ligation", "adapter_ligation"),
            ("Sample Cleanup 1", "sample_cleanup_1"),
            ("SPE Target Enrichment", "spe_target_enrichment"),
            ("Sample Cleanup 2", "sample_cleanup_2"),
            ("Universal PCR", "universal_pcr")
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
        self.lmgr = LayoutManager('QIAseq RNA Fusion XP Panels.lay')
        
        # Deck Resources - HSP Plates (only 2 available for general pipetting)
        self.HSP_Pipette = layout_item(self.lmgr, Plate96, 'HSP_Pipette')
        self.HSP_Pipette2 = layout_item(self.lmgr, Plate96, 'HSP_Pipette2')
        self.HSP_CPAC = layout_item(self.lmgr, Plate96, 'HSP_CPAC')
        self.HSP_Waste = layout_item(self.lmgr, Plate96, 'HSP_Waste')
        self.HSP_ODTC = layout_item(self.lmgr, Plate96, 'HSP_ODTC')
        self.HSP_ODTC_Lid = layout_item(self.lmgr, Lid, 'Ham_ComfortLid_ODTC')
        
        # MIDI Plates
        self.MIDI_Pipette = layout_item(self.lmgr, Plate96, 'MIDI_Pipette')
        self.MIDI_OnMagnet = layout_item(self.lmgr, Plate96, 'MIDI_OnMagnet')
        self.MIDI_Waste = layout_item(self.lmgr, Plate96, 'MIDI_Waste')
        
        # Small volume reagents (<2mL) - Use EppiCarrier32
        self.CAR_VIALS_SMALL = layout_item(self.lmgr, ReagentTrackedEppiCarrier32, 'CAR_VIALS_SMALL')
        
        # Assign reagent positions for small volumes
        self.FastSelect_position = self.CAR_VIALS_SMALL.assign_reagent_map('FastSelect', [11])
        self.RP_Primer_II_position = self.CAR_VIALS_SMALL.assign_reagent_map('RP_Primer_II', [1])
        self.First_Strand_Synthesis_Mix_position = self.CAR_VIALS_SMALL.assign_reagent_map('First_Strand_Mix', [2])
        self.Second_Strand_Synthesis_Mix_position = self.CAR_VIALS_SMALL.assign_reagent_map('Second_Strand_Mix', [3])
        self.ERAT_Mix_position = self.CAR_VIALS_SMALL.assign_reagent_map('ERAT_Mix', [4])
        self.Ligation_Mix_position = self.CAR_VIALS_SMALL.assign_reagent_map('Ligation_Mix', [6])
        self.UniversalPCR_position = self.CAR_VIALS_SMALL.assign_reagent_map('UniversalPCR', [10])
        
        # Bulk reagents (>5mL) - Use generic numbered reservoirs
        self.RGT_01 = layout_item(self.lmgr, ReagentTrackedReservoir60mL, 'RGT_01')  # Using existing name from layout
        self.RGT_02 = layout_item(self.lmgr, ReagentTrackedReservoir60mL, 'RGT_02')
        self.QIAseq_Beads_positions = self.RGT_01.assign_reagent_map('QIAseq_Beads', range(8))
        self.Nuclease_Free_Water_positions = self.RGT_01.assign_reagent_map('Nuclease_Free_Water', range(8))

        self.EthanolReservoir = layout_item(self.lmgr, ReagentTrackedBulkPlate, 'RGT_Ethanol')

        # SPE Master Mix - special carrier
        self.SPE_MasterMix_Container = layout_item(self.lmgr, ReagentTrackedReservoir60mL, 'RGT_03')
        self.SPE_MasterMix = self.SPE_MasterMix_Container.assign_reagent_map('SPE_MasterMix', range(8))
        
        
        # Stacked resources
        self.HSP_Stack = StackedResources.from_prefix(
            tracker_id="BioRadHardshell_Stack1",
            prefix="BioRadHardshell_Stack1",
            count=4,
            lmgr=self.lmgr,
            resource_type=Plate96
            )
        
        self.Lid_Stack = StackedResources.from_prefix(
            tracker_id="Ham_ComfortLid_ParkPos",
            prefix="Ham_ComfortLid_ParkPos",
            count=4,
            lmgr=self.lmgr,
            resource_type=Lid
            )
        
        self.MIDI_Stack = StackedResources.from_prefix(
            tracker_id="ABGENE_MIDI_Stack1",
            prefix="ABGENE_MIDI_Stack1",
            count=4,
            lmgr=self.lmgr,
            resource_type=Plate96
            )
        
        # Tracked tips
        self.tracked_tips_50uL = TrackedTips.from_prefix(
            tracker_id="TIP_50ulF_L",
            volume_capacity=50,
            prefix="TIP_50ulF_L",
            count=8,
            tip_type=Tip96,
            lmgr=self.lmgr)
        
        self.tracked_tips_300uL = TrackedTips.from_prefix(
            tracker_id="stf_l",
            volume_capacity=300,
            prefix="stf_l",
            count=8,
            tip_type=Tip96,
            lmgr=self.lmgr)

        self.tracked_tips_1000uL = TrackedTips.from_prefix(
            tracker_id="HTF_L",
            volume_capacity=1000,
            prefix="HTF_L",
            count=1,
            tip_type=Tip96,
            lmgr=self.lmgr)

        # Tip support tracker
        self.tip_support_resource = layout_item(self.lmgr, Tip96, 'TipSupport_0001')
        self.tip_support = TipSupportTracker(self.tip_support_resource)
        
        # Batch tracked objects for resource consumption
        self.tracked_reagent_vessels = [self.CAR_VIALS_SMALL, self.RGT_01, self.SPE_MasterMix_Container]
        self.tracked_tips = [self.tracked_tips_50uL, self.tracked_tips_300uL, self.tracked_tips_1000uL]
        self.stacked_resources = [self.HSP_Stack, self.MIDI_Stack]
    
    def _setup_volumes(self):
        """Setup volume parameters based on sample volume."""
        self.fast_select_volume = 50
        self.magbead_mix_volume = 1000
        self.post_shear_magbead_volume = self.sample_volume
        self.first_supernatant_removal_volume = self.sample_volume * 2
        self.supernatant_removal_volume = self.sample_volume + self.post_shear_magbead_volume
        self.m1_mix_volume = self.sample_volume * 1.6
        self.post_shear_etoh_wash_volume = 200
        self.post_shear_elution_buffer_volume = 30
        self.post_shear_elution_volume = 25.5
        
        # Calculate total reagent volumes needed (with 10% excess)
        total_samples_with_excess = self.num_samples * 1.1
        
        # Small volume reagents (<2mL total)
        self.fast_select_total = self.fast_select_volume * total_samples_with_excess / 1000  # mL
        self.rp_primer_total = 5 * total_samples_with_excess / 1000  # mL
        
        # Bulk reagents (>5mL total)
        self.ethanol_total = self.post_shear_etoh_wash_volume * 4 * total_samples_with_excess / 1000  # mL (2 washes x 2 cleanups)
        self.beads_total = self.post_shear_magbead_volume * 2 * total_samples_with_excess / 1000  # mL (2 cleanups)
        self.water_total = self.post_shear_elution_buffer_volume * 2 * total_samples_with_excess / 1000  # mL
    
    def _setup_hhs_devices(self):
        """Setup HHS (Hamilton Heater Shaker) devices."""
        self.HHS1_HSP = HHS(node=1, sequence="HHS1_HSP", lmgr=self.lmgr)
        self.HHS2_HSP = HHS(node=2, sequence="HHS2_HSP", lmgr=self.lmgr)
        self.HHS3_MIDI = HHS(node=3, sequence="HHS3_MIDI", lmgr=self.lmgr)
        self.HHS4_MIDI = HHS(node=4, sequence="HHS4_MIDI", lmgr=self.lmgr)
        self.HHS5_MIDI = HHS(node=5, sequence="HHS5_MIDI", lmgr=self.lmgr)

    def initialize(self):
        """Initialize the Hamilton system and all devices."""
        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())
            # Initialize HHS devices
            hhs_set_simulation(ham_int, 1 if self.device_simulation else 0)
            for node in range(1, 6):
                try:
                    hhs_create_usb_device(ham_int, node)
                    print(f"Initialized HHS node {node}")
                except Exception as e:
                    print(f"Warning: Could not initialize HHS node {node}: {e}")
            
            # Initialize CPAC for cold reagent storage
            initialize_cpac(ham_int, controller_id=1, simulating=self.device_simulation)
            set_temperature_target_cpac(ham_int, controller_id=1, device_id=1, target_temp=4.0)
            start_temperature_control_cpac(ham_int, controller_id=1, device_id=1)

            device_id = odtc_connect(ham_int, simulation_mode=self.device_simulation, local_ip='1.2.3.4', device_ip='5.6.7.8', device_port='COM4')  # Replace with actual IPs
            odtc_initialize(ham_int, device_id=device_id)

            print("Here's how many samples we are doing:", self.num_samples)

            
            self._initialized = True
            print("Protocol initialization completed successfully.")
    
    def first_strand_dna_synthesis(self):
        """Step 1: First Strand DNA Synthesis."""
        print("Starting First Strand DNA Synthesis...")

        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())
            
            # Transport HSP_Pipette2 to HSP_CPAC for reagent addition
            transport_resource(ham_int, self.HSP_Pipette2, self.HSP_CPAC, 
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Prepare sample positions
            HSP_CPAC_positions = [(self.HSP_CPAC, idx) for idx in range(self.num_samples)]
            volumes = [self.fast_select_volume] * self.num_samples
            
            # Add FastSelect reagent
            pip_transfer(ham_int, self.tracked_tips_50uL, self.FastSelect_position, 
                        HSP_CPAC_positions, volumes,
                        liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty',
                        aspiration_height=1, dispense_height=1)

            # Add RP Primer II
            pip_transfer(ham_int, self.tracked_tips_50uL, self.RP_Primer_II_position, 
                        HSP_CPAC_positions, volumes,
                        liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty',
                        aspiration_height=1, dispense_height=1)

            # Move plate from CPAC to HSP_Pipette2 for ODTC transfer
            transport_resource(ham_int, self.HSP_CPAC, self.HSP_Pipette2,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Move to thermal cycler
            odtc_open_door(ham_int, device_id=1)
            transport_resource(ham_int, self.HSP_Pipette2, self.HSP_ODTC,
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            
            # Add lid
            transport_resource(ham_int, self.Lid_Stack.fetch_next(), self.HSP_ODTC_Lid,
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.LID, iswap=True)

            odtc_close_door(ham_int, device_id=1)
            
            # Run thermal cycling protocol
            odtc_execute_protocol(ham_int, device_id=1, method_name='FirstStrandDNASynthesis', simulating=self.device_simulation)

            # Mix First Strand Synthesis Mix
            pip_mix(ham_int, tips=self.tracked_tips_50uL, positions_to_mix=self.First_Strand_Synthesis_Mix_position,
                   mix_volume=self.m1_mix_volume, mix_cycles=10,
                   liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty',
                   liquid_height=1)

            # Add First Strand Synthesis Mix to samples
            pip_transfer(ham_int, self.tracked_tips_50uL, self.First_Strand_Synthesis_Mix_position,
                        HSP_CPAC_positions, volumes,
                        liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty',
                        aspiration_height=1, dispense_height=1)


            
            
            # Wait for thermal cycler
            odtc_wait_for_idle(ham_int, device_id=1, simulating=self.device_simulation, check_interval=5)
            
            # Remove from thermal cycler
            odtc_open_door(ham_int, device_id=1)
            transport_resource(ham_int, self.HSP_ODTC_Lid, self.Lid_Stack.put_back(),
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.LID, iswap=True)
            transport_resource(ham_int, self.HSP_ODTC, self.HSP_Pipette2,
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            odtc_close_door(ham_int, device_id=1)


            print("First Strand DNA Synthesis completed.")
    
    def second_strand_dna_synthesis(self):
        """Step 2: Second Strand DNA Synthesis."""
        print("Starting Second Strand DNA Synthesis...")

        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())
            
            # Transport HSP_Pipette2 to HSP_CPAC
            transport_resource(ham_int, self.HSP_Pipette2, self.HSP_CPAC,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Mix Second Strand Synthesis Mix
            pip_mix(ham_int, tips=self.tracked_tips_50uL, positions_to_mix=self.Second_Strand_Synthesis_Mix_position,
                   mix_volume=self.m1_mix_volume, mix_cycles=10,
                   liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty',
                   liquid_height=1)
            
            # Add Second Strand Synthesis Mix
            HSP_CPAC_positions = [(self.HSP_CPAC, idx) for idx in range(self.num_samples)]
            volumes = [self.m1_mix_volume] * self.num_samples
            pip_transfer(ham_int, self.tracked_tips_300uL, self.Second_Strand_Synthesis_Mix_position,
                        HSP_CPAC_positions, volumes,
                        liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                        aspiration_height=1, dispense_height=1)

            # Move plate from CPAC to HSP_Pipette2
            transport_resource(ham_int, self.HSP_CPAC, self.HSP_Pipette2,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Move to thermal cycler
            odtc_open_door(ham_int, device_id=1)
            transport_resource(ham_int, self.HSP_Pipette2, self.HSP_ODTC,
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            
            # Add lid
            transport_resource(ham_int, self.Lid_Stack.fetch_next(), self.HSP_ODTC_Lid,
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.LID, iswap=True)
            
            odtc_close_door(ham_int, device_id=1)
            
            # Run thermal cycling protocol
            odtc_execute_protocol(ham_int, device_id=1, method_name='SecondStrandDNASynthesis', simulating=self.device_simulation)
            
            # Remove from thermal cycler
            odtc_open_door(ham_int, device_id=1)

            transport_resource(ham_int, self.HSP_ODTC_Lid, self.Lid_Stack.put_back(),
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.LID, iswap=True)

            transport_resource(ham_int, self.HSP_ODTC, self.HSP_Pipette2,
                                grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            
            odtc_close_door(ham_int, device_id=1)
            
            print("Second Strand DNA Synthesis completed.")
    
    def end_repair_a_tailing(self):
        """Step 3: End Repair and A-Tailing."""
        print("Starting End Repair and A-Tailing...")
        
        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())
            
            # Transport HSP_Pipette2 to HSP_CPAC
            transport_resource(ham_int, self.HSP_Pipette2, self.HSP_CPAC,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Mix ERAT Mix
            pip_mix(ham_int, self.tracked_tips_300uL, positions_to_mix=self.ERAT_Mix_position,
                   mix_volume=self.m1_mix_volume, mix_cycles=10,
                   liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                   liquid_height=1)
            
            # Add ERAT Mix
            HSP_CPAC_positions = [(self.HSP_CPAC, idx) for idx in range(self.num_samples)]
            volumes = [self.m1_mix_volume] * self.num_samples
            pip_transfer(ham_int, self.tracked_tips_300uL, self.ERAT_Mix_position,
                        HSP_CPAC_positions, volumes,
                        liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                        aspiration_height=1, dispense_height=1)

            # Move plate from CPAC to HSP_Pipette2
            transport_resource(ham_int, self.HSP_CPAC, self.HSP_Pipette2,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Move to thermal cycler
            odtc_open_door(ham_int, device_id=1)
            transport_resource(ham_int, self.HSP_Pipette2, self.HSP_ODTC,
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            
            # Add lid
            transport_resource(ham_int, self.Lid_Stack.fetch_next(), self.HSP_ODTC_Lid,
                             resource_type=GrippedResource.LID, core_gripper=True)
            
            odtc_close_door(ham_int, device_id=1)
            
            # Run thermal cycling protocol
            odtc_execute_protocol(ham_int, device_id=1, method_name='EndRepairATailing', simulating=self.device_simulation)

            # Remove from thermal cycler
            odtc_open_door(ham_int, device_id=1)
            transport_resource(ham_int, self.HSP_ODTC_Lid, self.Lid_Stack.put_back(),
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.LID, iswap=True)
            transport_resource(ham_int, self.HSP_ODTC, self.HSP_CPAC,
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            odtc_close_door(ham_int, device_id=1)
            
            print("End Repair and A-Tailing completed.")
    
    def adapter_ligation(self):
        """Step 4: Adapter Ligation."""
        print("Starting Adapter Ligation...")
        
        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())
            
            # Transport HSP_Pipette2 to HSP_CPAC
            transport_resource(ham_int, self.HSP_Pipette2, self.HSP_CPAC,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Prepare positions
            HSP_CPAC_positions = [(self.HSP_CPAC, idx) for idx in range(self.num_samples)]
            volumes = [self.sample_volume] * self.num_samples
            QIAseqIndexAdapter = self.HHS1_HSP.resource
            
            # Transfer QIAseq Index Adapters
            transfer_96(ham_int, tips=self.tracked_tips_50uL, tip_support=self.tip_support, source_plate=QIAseqIndexAdapter,
                        target_plate=self.HSP_CPAC, volume=self.sample_volume, num_samples=self.num_samples,
                        liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty',
                        aspiration_height=1)
            
            # Mix Ligation Mix
            pip_mix(ham_int, tips=self.tracked_tips_50uL, positions_to_mix=self.Ligation_Mix_position,
                   mix_volume=self.m1_mix_volume, mix_cycles=10,
                   liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty',
                   liquid_height=1)
            
            # Add Ligation Mix
            pip_transfer(ham_int, self.tracked_tips_50uL, self.Ligation_Mix_position,
                        HSP_CPAC_positions, volumes,
                        liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty',
                        aspiration_height=1, dispense_height=1)

            # Move plate from CPAC to HSP_Pipette2
            transport_resource(ham_int, self.HSP_CPAC, self.HSP_Pipette2,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Move to thermal cycler
            odtc_open_door(ham_int, device_id=1)
            transport_resource(ham_int, self.HSP_Pipette2, self.HSP_ODTC,
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            
            # Add lid
            transport_resource(ham_int, self.Lid_Stack.fetch_next(), self.HSP_ODTC_Lid,
                             resource_type=GrippedResource.LID, core_gripper=True)
            
            odtc_close_door(ham_int, device_id=1)
            
            # Run thermal cycling protocol
            odtc_execute_protocol(ham_int, device_id=1, method_name='AdapterLigation', simulating=self.device_simulation)

            # Remove from thermal cycler
            odtc_open_door(ham_int, device_id=1)
            transport_resource(ham_int, self.HSP_ODTC_Lid, self.Lid_Stack.put_back(),
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.LID, iswap=True)
            transport_resource(ham_int, self.HSP_ODTC, self.HSP_CPAC,
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            odtc_close_door(ham_int, device_id=1)
            
            print("Adapter Ligation completed.")
    
    def sample_cleanup_1(self):
        """Step 5: First Sample Cleanup with magnetic beads."""
        print("Starting Sample Cleanup 1...")

        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())
            
            # Get fresh MIDI plate from stack
            transport_resource(ham_int, self.MIDI_Stack.fetch_next(), self.HHS3_MIDI.resource,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Add QIAseq Beads to HHS position
            HHS3_MIDI_positions = [(self.HHS3_MIDI.resource, idx) for idx in range(self.num_samples)]
            volumes = [self.post_shear_magbead_volume] * self.num_samples
            pip_transfer(ham_int, self.tracked_tips_300uL, self.QIAseq_Beads_positions,
                        HHS3_MIDI_positions, volumes,
                        liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                        aspiration_height=1, dispense_height=1)

            transfer_96(ham_int, self.tracked_tips_300uL, tip_support=self.tip_support,num_samples=self.num_samples,
                        target_plate=self.HSP_CPAC, source_plate=self.HHS3_MIDI.resource, volume=self.sample_volume,
                        liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                        aspiration_height=1, dispense_height=1)
            
            # Shake plate
            hhs_start_shaker(ham_int, self.HHS3_MIDI.node, 1000)
            shake_timer = start_timer(10)
            shake_timer.wait(skip=self.device_simulation)
            hhs_stop_shaker(ham_int, self.HHS3_MIDI.node)
            
            # Transport to magnet
            transport_resource(ham_int, self.HHS3_MIDI.resource, self.MIDI_OnMagnet,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Let beads settle
            settle_timer = start_timer(60)
            settle_timer.wait(skip=self.device_simulation)
            

            # Ethanol wash (2x)
            for _ in range(2):
                ethanol_wash(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                            ethanol_plate=self.EthanolReservoir, magnet_plate=self.MIDI_OnMagnet,
                            waste_plate=self.MIDI_Waste, wash_volume=self.post_shear_etoh_wash_volume,
                            first_removal_volume=200, second_removal_volume=50,
                            liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

                        
            # Air dry
            dry_timer = start_timer(60)
            dry_timer.wait(skip=self.device_simulation)
            
            # Add Nuclease-Free Water
            MIDI_OnMagnet_positions = [(self.MIDI_OnMagnet, idx) for idx in range(self.num_samples)]
            pip_transfer(ham_int, self.tracked_tips_300uL, self.Nuclease_Free_Water_positions,
                        MIDI_OnMagnet_positions, [self.post_shear_elution_buffer_volume]*self.num_samples,
                        liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                        aspiration_height=1, dispense_height=1)

            # Transport to HHS3 to shake and incubate
            transport_resource(ham_int, self.MIDI_OnMagnet, self.HHS3_MIDI.resource,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Shake
            hhs_start_shaker(ham_int, self.HHS3_MIDI.node, 1000)
            shake_timer = start_timer(30)
            shake_timer.wait(skip=self.device_simulation)
            hhs_stop_shaker(ham_int, self.HHS3_MIDI.node)

            # Second round of beads - mix beads and add to HHS3 position
            HHS3_MIDI_positions = [(self.HHS3_MIDI.resource, idx) for idx in range(self.num_samples)]
            pip_mix(ham_int, tips=self.tracked_tips_300uL, positions_to_mix=self.QIAseq_Beads_positions,
                   mix_volume=self.post_shear_elution_volume, mix_cycles=10,
                   liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                   liquid_height=1)
            
            pip_transfer(ham_int, self.tracked_tips_300uL, self.QIAseq_Beads_positions,
                        HHS3_MIDI_positions, volumes,
                        liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                        aspiration_height=1, dispense_height=1)


            # Transport to magnet
            transport_resource(ham_int, self.HHS3_MIDI.resource, self.MIDI_OnMagnet,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Let beads settle
            settle_timer = start_timer(60)
            settle_timer.wait(skip=self.device_simulation)
            

            # Ethanol wash (2x)
            for _ in range(2):
                ethanol_wash(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                            ethanol_plate=self.EthanolReservoir, magnet_plate=self.MIDI_OnMagnet,
                            waste_plate=self.MIDI_Waste, wash_volume=self.post_shear_etoh_wash_volume,
                            first_removal_volume=200, second_removal_volume=50,
                            liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

                        
            # Air dry
            dry_timer = start_timer(60)
            dry_timer.wait(skip=self.device_simulation)

            # Final elution buffer
            pip_transfer(ham_int, self.tracked_tips_300uL, self.Nuclease_Free_Water_positions,
                        MIDI_OnMagnet_positions, [self.post_shear_elution_volume]*self.num_samples,
                        liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                        aspiration_height=1, dispense_height=1)


            # Bring HSP plate from stack to HSP_Pipette2
            transport_resource(ham_int, self.HSP_Stack.fetch_next(), self.HSP_Pipette2,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Move to HHS5 for mixing
            transport_resource(ham_int, self.MIDI_OnMagnet, self.HHS5_MIDI.resource,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Shake
            hhs_start_shaker(ham_int, self.HHS5_MIDI.node, 1000)
            shake_timer = start_timer(10)
            shake_timer.wait(skip=self.device_simulation)
            hhs_stop_shaker(ham_int, self.HHS5_MIDI.node)
            
            
            
            # Move to HHS5 for final mixing
            transport_resource(ham_int, self.MIDI_OnMagnet, self.HHS5_MIDI.resource,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Final shake
            hhs_start_shaker(ham_int, self.HHS5_MIDI.node, 1000)
            shake_timer = start_timer(10)
            shake_timer.wait(skip=self.device_simulation)
            hhs_stop_shaker(ham_int, self.HHS5_MIDI.node)
            
            # Back to magnet for elution
            transport_resource(ham_int, self.HHS5_MIDI.resource, self.MIDI_OnMagnet,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Get fresh HSP plate from stack
            transport_resource(ham_int, self.HSP_Stack.fetch_next(), self.HSP_Pipette2,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Transfer eluted samples to HSP plate
            HSP_Pipette2_positions = [(self.HSP_Pipette2, idx) for idx in range(self.num_samples)]
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                       self.MIDI_OnMagnet, self.HSP_Pipette2,
                       liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                       volume=self.post_shear_elution_volume,
                       aspiration_height=0.3, dispense_height=1)
            
            # Clean up - move used plates to waste
            transport_resource(ham_int, self.MIDI_OnMagnet, self.MIDI_Waste,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            transport_resource(ham_int, self.HSP_CPAC, self.HSP_Waste,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            print("Sample Cleanup 1 completed.")
    
    def spe_target_enrichment(self):
        """Step 6: SPE Target Enrichment."""
        print("Starting SPE Target Enrichment...")

        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())
            
            # Transport HSP_Pipette2 to HSP_CPAC
            transport_resource(ham_int, self.HSP_Pipette2, self.HSP_CPAC,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Mix SPE Master Mix
            pip_mix(ham_int, tips=self.tracked_tips_50uL, positions_to_mix=self.SPE_MasterMix,
                   mix_volume=self.m1_mix_volume, mix_cycles=10,
                   liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty',
                   liquid_height=1)
            
            # Add SPE Master Mix
            HSP_CPAC_positions = [(self.HSP_CPAC, idx) for idx in range(self.num_samples)]
            volumes = [self.m1_mix_volume] * self.num_samples
            pip_transfer(ham_int, self.tracked_tips_300uL, self.SPE_MasterMix,
                        HSP_CPAC_positions, volumes,
                        liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                        aspiration_height=1, dispense_height=1)

            # Move plate from CPAC to HSP_Pipette2
            transport_resource(ham_int, self.HSP_CPAC, self.HSP_Pipette2,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Move to thermal cycler
            odtc_open_door(ham_int, device_id=1)
            transport_resource(ham_int, self.HSP_Pipette2, self.HSP_ODTC,
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            
            # Add lid
            transport_resource(ham_int, self.Lid_Stack.fetch_next(), self.HSP_ODTC_Lid,
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.LID, iswap=True)

            odtc_close_door(ham_int, device_id=1)
            
            # Run thermal cycling protocol
            odtc_execute_protocol(ham_int, device_id=1, method_name='SPE_Target_Enrichment', simulating=self.device_simulation)

            # Remove from thermal cycler
            odtc_open_door(ham_int, device_id=1)
            transport_resource(ham_int, self.HSP_ODTC_Lid, self.Lid_Stack.put_back(),
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.LID, iswap=True)
            transport_resource(ham_int, self.HSP_ODTC, self.HSP_CPAC,
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            odtc_close_door(ham_int, device_id=1)
            
            print("SPE Target Enrichment completed.")
    
    def sample_cleanup_2(self):
        """Step 7: Second Sample Cleanup with magnetic beads."""
        print("Starting Sample Cleanup 2...")
        
        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())
            
            # Move HSP_Pipette2 to HSP_CPAC
            transport_resource(ham_int, self.HSP_Pipette2, self.HSP_CPAC,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Mix beads
            pip_mix(ham_int, tips=self.tracked_tips_300uL, positions_to_mix=self.QIAseq_Beads_positions,
                   mix_volume=self.magbead_mix_volume, mix_cycles=10,
                   liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                   liquid_height=1)

            # Dispense beads to MIDI_Pipette
            MIDI_Pipette_positions = [(self.MIDI_Pipette, idx) for idx in range(self.num_samples)]
            multi_dispense(ham_int, self.tracked_tips_300uL, self.QIAseq_Beads_positions,
                          MIDI_Pipette_positions, [self.post_shear_magbead_volume]*self.num_samples,
                          liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                          aspiration_height=1, dispense_height=1)

            # Add water
            multi_dispense(ham_int, self.tracked_tips_300uL, self.Nuclease_Free_Water_positions,
                          MIDI_Pipette_positions, [self.post_shear_elution_buffer_volume]*self.num_samples,
                          liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                          aspiration_height=1, dispense_height=1)

            # Transfer samples with 96 channel head
            HSP_Pipette2_positions = [(self.HSP_Pipette2, idx) for idx in range(self.num_samples)]
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                       self.HSP_Pipette2, self.MIDI_Pipette,
                       liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                       volume=self.post_shear_elution_volume,
                       aspiration_height=0.3, dispense_height=1)
            
            # Transport to HHS3 for mixing
            transport_resource(ham_int, self.MIDI_Pipette, self.HHS3_MIDI.resource,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Shake
            hhs_start_shaker(ham_int, self.HHS3_MIDI.node, 1000)
            shake_timer = start_timer(10)
            shake_timer.wait(skip=self.device_simulation)
            hhs_stop_shaker(ham_int, self.HHS3_MIDI.node)
            
            # Transport to magnet
            transport_resource(ham_int, self.HHS3_MIDI.resource, self.MIDI_OnMagnet,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Let beads settle
            settle_timer = start_timer(60)
            settle_timer.wait(skip=self.device_simulation)
            
            # Remove supernatant
            MIDI_OnMagnet_positions = [(self.MIDI_OnMagnet, idx) for idx in range(self.num_samples)]
            double_aspirate_supernatant_96(ham_int, self.tracked_tips_300uL, self.tip_support,
                                          self.num_samples, self.MIDI_OnMagnet, self.MIDI_Waste,
                                          self.first_supernatant_removal_volume, 30,
                                          liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                                          second_aspiration_height=0)
            
            # Ethanol wash (2x)
            for _ in range(2):
                ethanol_wash(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                            ethanol_plate=self.EthanolReservoir, magnet_plate=self.MIDI_OnMagnet,
                            waste_plate=self.MIDI_Waste, wash_volume=self.post_shear_etoh_wash_volume,
                            first_removal_volume=200, second_removal_volume=50,
                            liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')
            
            # Air dry
            dry_timer = start_timer(60)
            dry_timer.wait(skip=self.device_simulation)
            
            # Add elution buffer
            pip_transfer(ham_int, self.tracked_tips_300uL, self.Nuclease_Free_Water_positions,
                        MIDI_OnMagnet_positions, [self.post_shear_elution_volume]*self.num_samples,
                        liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                        aspiration_height=1, dispense_height=1)

            # Move to HHS5 for mixing
            transport_resource(ham_int, self.MIDI_OnMagnet, self.HHS5_MIDI.resource,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Shake
            hhs_start_shaker(ham_int, self.HHS5_MIDI.node, 1000)
            shake_timer = start_timer(10)
            shake_timer.wait(skip=self.device_simulation)
            hhs_stop_shaker(ham_int, self.HHS5_MIDI.node)
            
            # Back to magnet
            transport_resource(ham_int, self.HHS5_MIDI.resource, self.MIDI_OnMagnet,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Elute to HHS1_HSP
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                       self.MIDI_OnMagnet, self.HHS1_HSP.resource,
                       liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                       volume=self.post_shear_elution_volume,
                       aspiration_height=0.3, dispense_height=1)
            
            # Clean up
            transport_resource(ham_int, self.MIDI_OnMagnet, self.MIDI_Waste,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Move HHS1_HSP to thermal cycler for final enrichment
            odtc_open_door(ham_int, device_id=1)
            transport_resource(ham_int, self.HHS1_HSP.resource, self.HSP_ODTC,
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            
            # Add lid
            transport_resource(ham_int, self.Lid_Stack.fetch_next(), self.HSP_ODTC_Lid,
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.LID, iswap=True)
            
            odtc_close_door(ham_int, device_id=1)
            
            # Run final enrichment
            odtc_execute_protocol(ham_int, device_id=1, method_name='SecondStrandDNASynthesis', simulating=self.device_simulation)
            
            print("Sample Cleanup 2 completed.")
    
    def universal_pcr(self):
        """Step 8: Universal PCR."""
        print("Starting Universal PCR...")
        
        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())
            
            # Transport HSP_Pipette2 to HSP_CPAC
            transport_resource(ham_int, self.HSP_Pipette2, self.HSP_CPAC,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Add Universal PCR mix to adapter positions
            QIAseqIndexAdapter_positions = [(self.HHS1_HSP.resource, idx) for idx in range(self.num_samples)]
            print("Adding Universal PCR mix to adapter positions...")
            print(QIAseqIndexAdapter_positions)
            pip_transfer(ham_int, self.tracked_tips_300uL, self.UniversalPCR_position,
                        QIAseqIndexAdapter_positions, [self.m1_mix_volume]*self.num_samples,
                        liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                        aspiration_height=1, dispense_height=1)

            print("Universal PCR mix added.")
            # Mix
            mix_plate(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                     self.HHS1_HSP.resource, mixing_volume=self.m1_mix_volume, mix_cycles=10,
                     liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            # Get fresh HSP plate from stack
            transport_resource(ham_int, self.HSP_Stack.fetch_next(), self.HSP_CPAC,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Transfer Universal PCR mix to HSP_CPAC
            HSP_CPAC_positions = [(self.HSP_CPAC, idx) for idx in range(self.num_samples)]
            volumes = [self.m1_mix_volume] * self.num_samples
            pip_transfer(ham_int, self.tracked_tips_300uL, self.UniversalPCR_position,
                        HSP_CPAC_positions, volumes,
                        liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                        aspiration_height=1, dispense_height=1)

            # Move plate from CPAC to HSP_Pipette2
            transport_resource(ham_int, self.HSP_CPAC, self.HSP_Pipette2,
                             resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Move to thermal cycler
            odtc_open_door(ham_int, device_id=1)
            transport_resource(ham_int, self.HSP_Pipette2, self.HSP_ODTC,
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            
            # Add lid
            transport_resource(ham_int, self.Lid_Stack.fetch_next(), self.HSP_ODTC_Lid,
                             resource_type=GrippedResource.LID, grip_direction=GripDirection.RIGHT, iswap=True)
            
            odtc_close_door(ham_int, device_id=1)
            
            # Run Universal PCR protocol
            odtc_execute_protocol(ham_int, device_id=1, method_name='UniversalPCR', simulating=self.device_simulation)

            print("Universal PCR completed.")
    
    def run_complete_protocol(self):
        """Run the complete QIAseq RNA Fusion protocol."""
        print(f"Starting QIAseq RNA Fusion XP protocol for {self.num_samples} samples...")
        print(f"Sample volume: {self.sample_volume} µL")
        print(f"Simulation mode: {self.simulation}")
        
        try:
            # Run protocol steps
            self.initialize()
            self.first_strand_dna_synthesis()
            self.second_strand_dna_synthesis()
            self.end_repair_a_tailing()
            self.adapter_ligation()
            self.sample_cleanup_1()
            self.spe_target_enrichment()
            self.sample_cleanup_2()
            self.universal_pcr()
            
            print("QIAseq RNA Fusion XP protocol completed successfully!")
            
        except Exception as e:
            print(f"Protocol failed with error: {e}")
            raise


# Example usage
if __name__ == "__main__":
    # Create protocol instance
    protocol = QIAseqRNAFusionProtocol(
        num_samples=24,
        sample_volume=50,
        simulation=True,
        device_simulation=True
    )
    
    # Run complete protocol
    # protocol.run_complete_protocol()
    
    # Or run individual steps
    # protocol.initialize()
    # protocol.first_strand_dna_synthesis()
    # etc...
    
    # Generate consumables report
    protocol.run_protocol()
    generate_tadm_report()