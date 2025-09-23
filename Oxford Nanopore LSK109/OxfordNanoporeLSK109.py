# Oxford Nanopore LSK109 Protocol

from pyhamilton import (HamiltonInterface, LayoutManager, Plate96, Tip96, layout_item, normal_logging, start_timer)

from pyhamilton.pipetting import (pip_transfer, multi_dispense, double_aspirate_supernatant_96, ethanol_wash, transfer_96, pip_mix, mix_plate,
                                  shear_plate_96, transfer_96)
from pyhamilton.transport import (transport_resource, GripDirection, GrippedResource, GripperParams)
from pyhamilton.devices import (hhs_set_simulation, hhs_create_star_device, hhs_create_usb_device, hhs_set_temp_param,
                               hhs_start_temp_ctrl, hhs_stop_temp_ctrl, hhs_start_shaker, hhs_stop_shaker,
                                odtc_connect, odtc_initialize, odtc_open_door, odtc_close_door, odtc_execute_protocol,
                               initialize_cpac, set_temperature_target_cpac, start_temperature_control_cpac, get_temperature_cpac)
from pyhamilton.consumables import (tracked_volume_aspirate, ReagentTrackedReservoir60mL, ReagentTrackedPlate96, ReagentTrackedEppiCarrier32,
                                    ReagentTrackedPlate24, ReagentTrackedFalconCarrier24, generate_reagent_summary, ReagentTrackedBulkPlate)
from pyhamilton.ngs import Protocol, LoadingVis
from pyhamilton.resources import TrackedTips, StackedResources, Plate96, Tip96,Waste96, LayoutManager,TipSupportTracker, Lid, layout_item


import time
import os


class OxfordNanoporeLSK109Protocol(Protocol):
    """
    Oxford Nanopore LSK109 library preparation protocol class for Hamilton liquid handler.
    
    This class encapsulates the complete workflow for Oxford Nanopore LSK109 library preparation
    including DNA end prep, cleanup, adapter ligation, and final cleanup.
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
            ("cDNA End Prep", "cdna_end_prep"),
            ("End Prep Cleanup", "end_prep_cleanup"),
            ("Adapter Ligation", "adapter_ligation"),
            ("Adapter Ligation Cleanup", "adapter_ligation_cleanup")
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
        self._initialized = False
    
    def _setup_deck_layout(self):
        """Setup deck layout and resources."""
        # Deck layout import
        self.lmgr = LayoutManager('LSK109_MPH_OxfordNanopore_v0.9.2.lay')

        # HSP Deck Resources - only two available for general pipetting
        self.HSP_Pipette = layout_item(self.lmgr, Plate96, 'HSP_Pipette')
        self.HSP_Pipette2 = layout_item(self.lmgr, Plate96, 'HSP_Pipette2')  
        self.HSP_Waste = layout_item(self.lmgr, Plate96, 'HSP_Waste')
        self.HSP_ODTC = layout_item(self.lmgr, Plate96, 'HSP_ODTC')
        self.HSP_ODTC_Lid = layout_item(self.lmgr, Lid, 'Ham_ComfortLid_ODTC')
        
        # MIDI Deck Resources
        self.MIDI_Pipette = layout_item(self.lmgr, Plate96, 'MIDI_Pipette')  
        self.MIDI_OnMagnet = layout_item(self.lmgr, Plate96, 'MIDI_OnMagnet') 
        self.MIDI_Waste = layout_item(self.lmgr, Plate96, 'MIDI_Waste')
        self.HHS5_MIDI = layout_item(self.lmgr, Plate96, 'HHS5_MIDI')

        # Cold-sensitive reagent containers (enzymes) - use CPAC with temperature control
        self.CPAC_Reagents = layout_item(self.lmgr, ReagentTrackedPlate96, 'CPAC_Adapter')
        
        # Room temperature small volume reagents
        self.CAR_VIALS_SMALL = layout_item(self.lmgr, ReagentTrackedEppiCarrier32, 'CAR_VIALS_SMALL')

        # Bulk reagents - use generic numbered reservoirs
        self.RGT_01 = layout_item(self.lmgr, ReagentTrackedReservoir60mL, 'RGT_01')  # AMPure beads
        self.RGT_02 = layout_item(self.lmgr, ReagentTrackedReservoir60mL, 'RGT_02')  # Fragment Buffer
        self.RGT_03 = layout_item(self.lmgr, ReagentTrackedReservoir60mL, 'RGT_03')  # Nuclease Free Water
        self.EthanolReservoir = layout_item(self.lmgr, ReagentTrackedBulkPlate, 'Ethanol_Reservoir')  # 80% Ethanol

        # Reagent position assignments
        # Cold-sensitive reagents (CPAC)
        self.endprep_mix_positions = self.CPAC_Reagents.assign_reagent_map('EndPrepMix', [0])
        self.adapter_ligation_mix_positions = self.CPAC_Reagents.assign_reagent_map('AdapterLigationMix', [1])

        # Small volume room temperature reagents
        self.adapter_mix_positions = self.CAR_VIALS_SMALL.assign_reagent_map('AdapterMix', [0])
        self.elution_buffer_positions = self.CAR_VIALS_SMALL.assign_reagent_map('ElutionBuffer', [1])

        # Bulk reagents (8-channel access)
        self.ampure_bead_positions = self.RGT_01.assign_reagent_map('AMPureBeads', list(range(8)))
        self.fragment_buffer_positions = self.RGT_02.assign_reagent_map('FragmentBuffer', list(range(8)))
        self.nuclease_free_water_positions = self.RGT_03.assign_reagent_map('NucleaseFreeWater', list(range(8)))

        # Stacked resources
        self.HSP_Stack = StackedResources.from_prefix(
            tracker_id="BioRadHardshell_Stack1",
            prefix="BioRadHardshell_Stack1",
            lmgr=self.lmgr,
            resource_type=Plate96,
            count=4)

        self.MIDI_Stack = StackedResources.from_prefix(
            tracker_id="ABGENE_MIDI_Stack1",
            prefix="ABGENE_MIDI_Stack1",
            lmgr=self.lmgr,
            resource_type=Plate96,
            count=3)
        
        self.Lid_Stack = StackedResources.from_prefix(
            tracker_id="Ham_ComfortLid_Stack_ParkPos",
            prefix="Ham_ComfortLid_Stack_ParkPos",
            lmgr=self.lmgr,
            resource_type=Lid,
            count=4)

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


        self.tip_support_resource = layout_item(self.lmgr, Tip96, 'TipSupport_0001')
        self.tip_support = TipSupportTracker(self.tip_support_resource)

        # Batch together tracked objects for resource consumption estimates
        self.tracked_reagent_vessels = [
            self.CPAC_Reagents, self.CAR_VIALS_SMALL, self.RGT_01, self.RGT_02, self.RGT_03
        ]

        self.tracked_tips = [self.tracked_tips_50uL, self.tracked_tips_300uL]
        self.stacked_resources = [self.HSP_Stack, self.MIDI_Stack]

    def _setup_volumes(self):
        """Setup volume parameters based on sample volume and calculate total reagent needs."""
        # Per-sample volumes
        self.endprep_mix_volume = 7.5  # NEBNext Ultra II End Prep reaction buffer + enzyme per sample
        self.adapter_ligation_mix_volume = 20  # Ligation Master Mix per sample
        self.adapter_mix_volume = 5   # Adapter Mix per sample
        self.magbead_volume = self.sample_volume  # 1:1 ratio with sample
        self.ethanol_wash_volume = 200  # Per wash
        self.fragment_buffer_volume = 60  # Per wash
        self.elution_volume = 15  # Final elution volume
        
        # Supernatant removal volumes
        self.first_supernatant_removal_volume = self.sample_volume + self.magbead_volume
        self.supernatant_removal_volume = self.sample_volume + self.magbead_volume + 10  # Extra for safety
        
        # Calculate total reagent volumes needed (with 10% excess)
        excess_factor = 1.1
        self.total_endprep_mix_needed = self.endprep_mix_volume * self.num_samples * excess_factor
        self.total_adapter_ligation_mix_needed = self.adapter_ligation_mix_volume * self.num_samples * excess_factor
        self.total_adapter_mix_needed = self.adapter_mix_volume * self.num_samples * excess_factor
        self.total_magbead_needed = self.magbead_volume * self.num_samples * 3 * excess_factor  # 3 cleanup steps
        self.total_ethanol_needed = self.ethanol_wash_volume * self.num_samples * 3 * excess_factor  # 3 washes
        self.total_fragment_buffer_needed = self.fragment_buffer_volume * self.num_samples * 2 * excess_factor  # 2 washes
        self.total_elution_needed = self.elution_volume * self.num_samples * excess_factor

    def initialize(self):
        """Initialize the Hamilton system and all devices including CPAC temperature control."""
        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())
            
            hhs_set_simulation(ham_int, 1 if self.device_simulation else 0)
            for node in range(1, 6):
                try:
                    hhs_create_usb_device(ham_int, node)
                    print(f"Initialized HHS node {node}")
                except Exception as e:
                    print(f"Warning: Could not initialize HHS node {node}: {e}")
            
            # Initialize CPAC for cold reagents
            initialize_cpac(ham_int, controller_id=1, simulating = self.device_simulation)
            set_temperature_target_cpac(ham_int, 4.0, controller_id=1, device_id=1)  # 4°C for enzyme storage
            start_temperature_control_cpac(ham_int, controller_id=1, device_id=1)
            
            device_id = odtc_connect(ham_int, simulation_mode=self.device_simulation, local_ip='1.2.3.4', device_ip='5.6.7.8', device_port='COM4')  # Replace with actual IPs
            odtc_initialize(ham_int, device_id=device_id)

            
            self._initialized = True
            print("LSK109 Protocol initialization completed successfully.")
            print(f"Total reagent volumes calculated:")
            print(f"  End Prep Mix: {self.total_endprep_mix_needed:.1f} µL")
            print(f"  Adapter Ligation Mix: {self.total_adapter_ligation_mix_needed:.1f} µL")
            print(f"  Adapter Mix: {self.total_adapter_mix_needed:.1f} µL")
            print(f"  AMPure beads: {self.total_magbead_needed:.1f} µL")
            print(f"  80% Ethanol: {self.total_ethanol_needed:.1f} µL")
            print(f"  Fragment Buffer: {self.total_fragment_buffer_needed:.1f} µL")

    def cdna_end_prep(self):
        """Step 1: cDNA End Preparation."""
        if not self._initialized:
            print("Warning: Protocol not initialized. Call initialize() first.")
        
        print("Starting cDNA End Preparation step...")
        
        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())

            # Transfer EndPrep Master Mix to HSP plate (samples already loaded)
            hsp_positions = [(self.HSP_Pipette, i) for i in range(self.num_samples)]
            pip_transfer(ham_int, tips=self.tracked_tips_300uL, source_positions=self.endprep_mix_positions, 
                        dispense_positions=hsp_positions, volumes=[self.endprep_mix_volume] * self.num_samples,
                        liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                        aspiration_height=1, dispense_height=1)

            # Mix plate
            print("Tip support at first step")
            print(self.tip_support.source_rack)
            mix_plate(ham_int, self.tracked_tips_50uL, self.tip_support, self.num_samples, 
                     plate = self.HSP_Pipette, mixing_volume=25, mix_cycles=5,
                     liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')

            odtc_open_door(ham_int, device_id=1)

            # Move HSP plate to thermal cycler and run protocol
            transport_resource(ham_int, self.HSP_Pipette, self.HSP_ODTC, 
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)

            transport_resource(ham_int, self.Lid_Stack.fetch_next(), self.HSP_ODTC_Lid, 
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.LID, iswap=True)

            # User loads protocol on ODTC
            odtc_close_door(ham_int, device_id=1)

            odtc_execute_protocol(ham_int, device_id=1, method_name='LSK109_EndPrep', simulating=self.device_simulation)

            odtc_open_door(ham_int, device_id=1)

            # Remove lid and place back in stack
            transport_resource(ham_int, self.HSP_ODTC_Lid, self.Lid_Stack.put_back(),
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.LID, iswap=True)

            # Move plate back to HSP_Pipette2 position
            transport_resource(ham_int, self.HSP_ODTC, self.HSP_Pipette2, 
                             grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)

            print("cDNA End Preparation completed.")

    def end_prep_cleanup(self):
        """Step 2: End Prep Cleanup with magnetic beads."""
        print("Starting End Prep Cleanup step...")

        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())

            # Get fresh MIDI plate from stack
            transport_resource(ham_int, self.MIDI_Stack.fetch_next(), self.MIDI_Pipette, 
                             core_gripper=True, resource_type=GrippedResource.MIDI)

            # Transfer AMPure beads to MIDI plate
            midi_positions = [(self.MIDI_Pipette, i) for i in range(self.num_samples)]
            pip_transfer(ham_int, self.tracked_tips_50uL, self.ampure_bead_positions, midi_positions,
                        volumes=[self.magbead_volume] * self.num_samples,
                        liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty',
                        aspiration_height=1, dispense_height=1)

            # Transfer samples from HSP to MIDI plate
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                       self.HSP_Pipette2, self.MIDI_Pipette,
                       volume=self.sample_volume + self.endprep_mix_volume,
                       liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                       aspiration_height=1, dispense_height=1)

            # Move HSP plate to waste to free up position
            transport_resource(ham_int, self.HSP_Pipette2, self.HSP_Waste, 
                             core_gripper=True, resource_type=GrippedResource.PCR)

            # Mix beads and samples, then incubate
            mix_plate(ham_int, self.tracked_tips_50uL, self.tip_support, self.num_samples,
                     self.MIDI_Pipette, mix_cycles=5, mixing_volume=50,
                     liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')

            # Incubate at room temperature for 5 minutes
            incubation_timer = start_timer(300)
            incubation_timer.wait(skip=self.device_simulation)

            # Move to magnetic module
            transport_resource(ham_int, self.MIDI_Pipette, self.MIDI_OnMagnet, 
                             core_gripper=True, resource_type=GrippedResource.PCR)

            # Wait for beads to settle
            settle_timer = start_timer(120)
            settle_timer.wait(skip=self.device_simulation)

            # Remove supernatant
            double_aspirate_supernatant_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                                          self.MIDI_OnMagnet, self.MIDI_Waste,
                                          first_volume=self.first_supernatant_removal_volume,
                                          second_volume=20,
                                          liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                                          first_aspiration_height=0.5)

            # 2x Ethanol washes
            for wash_num in range(2):
                ethanol_wash(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                            ethanol_plate=self.EthanolReservoir, magnet_plate=self.MIDI_OnMagnet, 
                            waste_plate=self.MIDI_Waste, wash_volume=self.ethanol_wash_volume,
                            first_removal_volume=self.ethanol_wash_volume + 20,
                            second_removal_volume=30,
                            liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            # Air dry beads for 2 minutes
            dry_timer = start_timer(120)
            dry_timer.wait(skip=self.device_simulation)

            # Get fresh HSP plate for elution
            transport_resource(ham_int, self.HSP_Stack.fetch_next(), self.HSP_Pipette, 
                             core_gripper=True, resource_type=GrippedResource.PCR)

            # Elute with Nuclease Free Water
            MIDI_OnMagnet_positions = [(self.MIDI_OnMagnet, i) for i in range(self.num_samples)]
            pip_transfer(ham_int, tips=self.tracked_tips_50uL, 
                       source_positions=self.nuclease_free_water_positions, dispense_positions=MIDI_OnMagnet_positions,
                       volumes=[self.elution_volume] * self.num_samples,
                       liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty',
                       aspiration_height=1, dispense_height=1)

            # Mix and incubate
            mix_plate(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                     self.MIDI_OnMagnet, mixing_volume=12,mix_cycles=5,
                     liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            mix_timer = start_timer(120)
            mix_timer.wait(skip=self.device_simulation)

            # Transfer eluted DNA to HSP plate
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                       self.MIDI_OnMagnet, self.HSP_Pipette,
                       volume=self.elution_volume - 2,  # Leave 2µL to avoid beads
                       liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                       aspiration_height=3, dispense_height=1)

            # Move MIDI plate to waste
            transport_resource(ham_int, self.MIDI_OnMagnet, self.MIDI_Waste,
                             core_gripper=True, resource_type=GrippedResource.MIDI)

            print("End Prep Cleanup completed.")

    def adapter_ligation(self):
        """Step 3: Adapter Ligation."""
        print("Starting Adapter Ligation step...")

        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())

            # Add Adapter Ligation Master Mix to samples
            hsp_positions = [(self.HSP_Pipette, i) for i in range(self.num_samples)]
            pip_transfer(ham_int, self.tracked_tips_300uL, self.adapter_ligation_mix_positions, hsp_positions,
                        volumes=[self.adapter_ligation_mix_volume] * self.num_samples,
                        liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                        aspiration_height=1, dispense_height=1)

            # Add Adapter Mix
            pip_transfer(ham_int, self.tracked_tips_300uL, self.adapter_mix_positions, hsp_positions,
                        volumes=[self.adapter_mix_volume] * self.num_samples,
                        liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                        aspiration_height=1, dispense_height=1)

            # Mix plate
            mix_plate(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                     self.HSP_Pipette, mixing_volume=30,mix_cycles=5,
                     liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            # Incubate at room temperature for 10 minutes
            ligation_timer = start_timer(600)
            ligation_timer.wait(skip=self.device_simulation)

            print("Adapter Ligation completed.")

    def adapter_ligation_cleanup(self):
        """Step 4: Adapter Ligation Cleanup."""
        print("Starting Adapter Ligation Cleanup step...")

        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())

            # Get fresh MIDI plate
            transport_resource(ham_int, self.MIDI_Stack.fetch_next(), self.MIDI_Pipette,
                             core_gripper=True, resource_type=GrippedResource.PCR)

            # Add AMPure beads to MIDI plate
            midi_positions = [(self.MIDI_Pipette, i) for i in range(self.num_samples)]
            pip_transfer(ham_int, self.tracked_tips_300uL, self.ampure_bead_positions, midi_positions,
                        volumes=[self.magbead_volume] * self.num_samples,
                        liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                        aspiration_height=1, dispense_height=1)

            # Transfer ligation reaction to MIDI plate
            total_ligation_volume = self.elution_volume - 2 + self.adapter_ligation_mix_volume + self.adapter_mix_volume
            transfer_96(ham_int, tips=self.tracked_tips_300uL, tip_support=self.tip_support, num_samples=self.num_samples,
                       source_plate=self.HSP_Pipette, target_plate=self.MIDI_Pipette,
                       volume=total_ligation_volume,
                       liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                       aspiration_height=1, dispense_height=1)

            # Move HSP plate to waste
            transport_resource(ham_int, self.HSP_Pipette, self.HSP_Waste, 
                             core_gripper=True, resource_type=GrippedResource.PCR)

            # Mix and incubate
            mix_plate(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                     self.MIDI_Pipette, mixing_volume=50, mix_cycles=5,
                     liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            incubation_timer = start_timer(300)
            incubation_timer.wait(skip=self.device_simulation)

            # Move to magnet
            transport_resource(ham_int, self.MIDI_Pipette, self.MIDI_OnMagnet, 
                             core_gripper=True, resource_type=GrippedResource.PCR)

            settle_timer = start_timer(120)
            settle_timer.wait(skip=self.device_simulation)

            # Remove supernatant
            double_aspirate_supernatant_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                                          self.MIDI_OnMagnet, self.MIDI_Waste,
                                          first_volume=self.supernatant_removal_volume,
                                          second_volume=30,
                                          liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                                          first_aspiration_height=0.5)

            # 2x Fragment buffer washes
            for wash_cycle in range(2):
                # Add fragment buffer
                transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                           self.EthanolReservoir, self.MIDI_OnMagnet,
                           volume=self.fragment_buffer_volume,
                           liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                           aspiration_height=1, dispense_height=1)

                # Mix
                mix_plate(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                          self.MIDI_OnMagnet, mixing_volume=40,mix_cycles=5,
                          liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

                # Incubate
                wash_timer = start_timer(120)
                wash_timer.wait(skip=self.device_simulation)

                # Remove supernatant
                double_aspirate_supernatant_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                                              self.MIDI_OnMagnet, self.MIDI_Waste,
                                              first_volume=self.fragment_buffer_volume + 20,
                                              second_volume=30,
                                              liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                                              first_aspiration_height=0.5)

            # Air dry beads
            dry_timer = start_timer(120)
            dry_timer.wait(skip=self.device_simulation)

            # Get fresh HSP plate for final elution
            transport_resource(ham_int, self.HSP_Stack.fetch_next(), self.HSP_Pipette2, 
                             core_gripper=True, resource_type=GrippedResource.PCR)

            # Final elution with Elution Buffer
            MIDI_OnMagnet_positions = [(self.MIDI_OnMagnet, i) for i in range(self.num_samples)]
            pip_transfer(ham_int, self.tracked_tips_300uL,
                       self.elution_buffer_positions, MIDI_OnMagnet_positions,
                       volumes=[self.elution_volume] * self.num_samples,
                       liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                       aspiration_height=1, dispense_height=1)

            # Mix and incubate
            mix_plate(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                     self.MIDI_OnMagnet, mixing_volume=12, mix_cycles=5,
                     liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            final_mix_timer = start_timer(120)
            final_mix_timer.wait(skip=self.device_simulation)

            # Collect final eluted DNA
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                       self.MIDI_OnMagnet, self.HSP_Pipette2,
                       volume=self.elution_volume - 2,  # Leave 2µL to avoid beads
                       liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                       aspiration_height=0.5, dispense_height=1)

            # Move MIDI plate to waste
            transport_resource(ham_int, self.MIDI_OnMagnet, self.MIDI_Waste, 
                             core_gripper=True, resource_type=GrippedResource.PCR)

            print("Adapter Ligation Cleanup completed.")
            print(f"Final LSK109 libraries ready in {self.HSP_Pipette2.layout_name()}")

    def run_complete_protocol(self):
        """Run the complete Oxford Nanopore LSK109 protocol."""
        print(f"Starting Oxford Nanopore LSK109 protocol for {self.num_samples} samples...")
        print(f"Sample volume: {self.sample_volume} µL")
        print(f"Simulation mode: {self.simulation}")
        
        try:
            # Run protocol steps
            self.initialize()
            self.cdna_end_prep()
            self.end_prep_cleanup()
            self.adapter_ligation()
            self.adapter_ligation_cleanup()
            
            print("Oxford Nanopore LSK109 protocol completed successfully!")
            print("Libraries are ready for loading onto nanopore flow cells.")
            
        except Exception as e:
            print(f"Protocol failed with error: {e}")
            raise


# Example usage
if __name__ == "__main__":
    # Create protocol instance
    protocol = OxfordNanoporeLSK109Protocol(
        num_samples=96,
        sample_volume=50,
        simulation=True,
        device_simulation=True
    )

    # Run complete protocol
    # protocol.run_complete_protocol()

    # Or run individual steps
    # protocol.initialize()
    # protocol.cdna_end_prep()
    # protocol.end_prep_cleanup()
    # protocol.adapter_ligation()
    # protocol.adapter_ligation_cleanup()

    # Calculate consumables
    protocol.run_protocol()