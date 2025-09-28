from pyhamilton import (HamiltonInterface, LayoutManager, Plate96, Tip96, start_timer,
                        normal_logging, layout_item)

from pyhamilton.devices import (initialize_cpac, set_temperature_target_cpac, start_temperature_control_cpac,
                                 get_temperature_cpac, odtc_initialize, odtc_connect, odtc_open_door, odtc_close_door, odtc_execute_protocol,
                                 odtc_wait_for_idle, odtc_get_status)

from pyhamilton.consumables import (tracked_volume_aspirate, ReagentTrackedReservoir60mL, ReagentTrackedPlate96, 
                                    ReagentTrackedFalconCarrier24, generate_reagent_summary, ReagentTrackedBulkPlate, ReagentTrackedEppiCarrier32)

from pyhamilton.ngs import Protocol, LoadingVis, generate_tadm_report

from pyhamilton.resources import TrackedTips, StackedResources, Plate96, Tip96,Waste96, LayoutManager,TipSupportTracker, Lid, layout_item

from pyhamilton.pipetting import (pip_transfer, multi_dispense, double_aspirate_supernatant_96, ethanol_wash, transfer_96, pip_mix, mix_plate,
                                  shear_plate_96, transfer_96)
from pyhamilton.transport import (transport_resource, GripDirection, GrippedResource, GripperParams)
from pyhamilton.odtc.odtc_protocol import ThermalCyclerProtocol
import time
import os


class TenXGEXLibraryPrepProtocol(Protocol):
    """
    10X Genomics GEX Library Preparation protocol class for Hamilton liquid handler.
    
    This class encapsulates the complete workflow for 10X GEX library preparation
    including fragmentation, end repair, A-tailing, adapter ligation, cleanup,
    sample indexing PCR, and final size selection.
    """
    
    def __init__(self, num_samples=8, sample_volume=10, pcr_cycles=12, simulation=False, device_simulation=False):
        """
        Initialize the protocol with sample parameters.
        
        Args:
            num_samples (int): Number of samples to process (default: 8)
            sample_volume (int): Volume of each sample in µL (default: 10)
            pcr_cycles (int): Number of PCR cycles based on cDNA input (default: 12)
            simulation (bool): Whether to run in simulation mode (default: False)
        """

        self.available_steps = [
            ("Initialize System", "initialize"),
            ("Fragmentation, End Repair & A-tailing", "fragmentation_end_repair_atailing"),
            ("Post Fragmentation SPRIselect", "post_fragmentation_spriselect"),
            ("Adapter Ligation", "adapter_ligation"),
            ("Post Ligation Cleanup", "post_ligation_cleanup"),
            ("Sample Index PCR", "sample_index_pcr"),
            ("Final Size Selection", "final_size_selection")
        ]

        super().__init__()
        self.num_samples = num_samples
        self.sample_volume = sample_volume
        self.pcr_cycles = pcr_cycles
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
        self.lmgr = LayoutManager('10X_GEX_LibraryPrep_Deck.lay')

        # PCR plates for thermal cycling steps
        self.HSP_Plate = layout_item(self.lmgr, Plate96, 'HSP_Pipette')
        self.HSP_Plate2 = layout_item(self.lmgr, Plate96, 'HSP_Pipette2')
        self.HSP_Waste = layout_item(self.lmgr, Plate96, 'HSP_Waste')
        self.HHS1_HSP = layout_item(self.lmgr, Plate96, 'HHS1_HSP')
        self.HSP_Magnet = layout_item(self.lmgr, Plate96, 'HSP_OnMagnet')
        self.HSP_ODTC = layout_item(self.lmgr, Plate96, 'HSP_ODTC')
        self.HSP_ODTC_Lid = layout_item(self.lmgr, Lid, 'Ham_ComfortLid_ODTC')

        # MIDI plates for magnetic bead operations
        self.MIDI_Plate = layout_item(self.lmgr, Plate96, 'MIDI_Pipette')
        self.MIDI_Waste = layout_item(self.lmgr, Plate96, 'MIDI_Waste')
        self.MIDI_CPAC = layout_item(self.lmgr, Plate96, 'MIDI_CPAC')
        self.MIDI_OnMagnet = layout_item(self.lmgr, Plate96, 'MIDI_OnMagnet')
        self.HHS5_MIDI = layout_item(self.lmgr, Plate96, 'HHS5_MIDI')


        # Reagent containers
        self.CAR_VIALS_SMALL = layout_item(self.lmgr, ReagentTrackedEppiCarrier32, 'CAR_VIALS_SMALL')
        self.CPAC_Reagents = layout_item(self.lmgr, ReagentTrackedPlate96, 'CPAC_HSP_0001')
        self.RGT_01 = layout_item(self.lmgr, ReagentTrackedReservoir60mL, 'rgt_cont_60ml_BC_A00_0001')
        self.RGT_02 = layout_item(self.lmgr, ReagentTrackedReservoir60mL, 'rgt_cont_60ml_BC_A00_0002')
        self.EthanolReservoir = layout_item(self.lmgr, ReagentTrackedBulkPlate, 'RGT_Ethanol')

        # Reagent mapping and positions
        # Room temperature small volume reagents
        self.fragmentation_buffer_positions = self.CAR_VIALS_SMALL.assign_reagent_map('FragmentationBuffer', [0])
        self.ligation_mix_positions = self.CAR_VIALS_SMALL.assign_reagent_map('LigationMix', [1])
        self.buffer_eb_positions = self.CAR_VIALS_SMALL.assign_reagent_map('BufferEB', [2])

        # Cold-sensitive reagents in CPAC
        self.fragmentation_enzyme_positions = self.CPAC_Reagents.assign_reagent_map('FragmentationEnzyme', [0])
        self.dna_ligase_positions = self.CPAC_Reagents.assign_reagent_map('DNALigase', [1])
        self.library_amp_mix_positions = self.CPAC_Reagents.assign_reagent_map('LibraryAmpMix', [2])

        # Bulk reagents
        self.spriselect_positions = self.RGT_02.assign_reagent_map('SPRIselect', range(8))

        # Index positions
        self.index_positions = [(self.HHS1_HSP, i) for i in range(self.num_samples)]

        # Liquid waste
        self.Liquid_Waste = layout_item(self.lmgr, Waste96, 'core96externalwaste_0001')

        # Stacked resources
        self.HSP_Stack = StackedResources.from_prefix(
            tracker_id="BioRadHardShell_Stack4",
            prefix="BioRadHardShell_Stack4",
            lmgr=self.lmgr,
            resource_type=Plate96,
            count=5)

        self.Lid_Stack = StackedResources.from_prefix(
            tracker_id="Ham_ComfortLid_Stack", 
            prefix="Ham_ComfortLid_Stack",
            lmgr=self.lmgr,
            resource_type=Lid,
            count=3)
        
        self.MIDI_Stack = StackedResources.from_prefix(
            tracker_id="AbgeneMIDI_Stack1",
            prefix="AbgeneMIDI_Stack1",
            lmgr=self.lmgr,
            resource_type=Plate96,
            count=3)

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

        # Tip support
        self.tip_support_resource = layout_item(self.lmgr, Tip96, 'TipSupport_0001')
        self.tip_support = TipSupportTracker(self.tip_support_resource)

        # Batch together tracked objects for resource consumption estimates
        self.tracked_reagent_vessels = [
            self.CAR_VIALS_SMALL, self.CPAC_Reagents, self.RGT_01, self.RGT_02, self.EthanolReservoir
        ]

        self.tracked_tips = [self.tracked_tips_50uL, self.tracked_tips_300uL]
        self.stacked_resources = [self.HSP_Stack, self.Lid_Stack, self.MIDI_Stack]

    def _setup_volumes(self):
        """Setup volume parameters based on sample volume."""
        self.fragmentation_mix_volume = 40  # µL per sample
        self.ligation_mix_volume = 50       # µL per sample
        self.amp_mix_volume = 50            # µL per sample
        self.index_volume = 20              # µL per sample
        
        # Calculated volumes with excess
        self.excess_factor = 1.1

    def initialize(self):
        """Initialize Hamilton system and CPAC for cold reagents."""
        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())
            
            # Initialize CPAC for cold reagent storage
            initialize_cpac(ham_int, controller_id=1, simulating=self.simulation)
            set_temperature_target_cpac(ham_int, controller_id=1, device_id=1, target_temp=4)
            start_temperature_control_cpac(ham_int, controller_id=1, device_id=1)
            
            # Verify CPAC temperature
            temp = get_temperature_cpac(ham_int, controller_id=1, device_id=1)
            print(f"CPAC temperature stabilized at: {temp}°C")

            device_id = odtc_connect(ham_int, simulation_mode=self.device_simulation, local_ip='1.2.3.4', device_ip='5.6.7.8', device_port='COM4')  # Replace with actual IPs
            odtc_initialize(ham_int, device_id=device_id)

            
            self._initialized = True
            print("10X GEX Library Prep Protocol initialization completed successfully.")

    def _prepare_fragmentation_mix(self, ham_int):
        """Prepare fragmentation mix with 10% excess."""
        total_reactions = int(self.num_samples * self.excess_factor)
        
        buffer_eb_vol = 20 * total_reactions
        frag_buffer_vol = 5 * total_reactions
        frag_enzyme_vol = 10 * total_reactions
        
        # Mix fragmentation buffer first (room temp reagent)
        pip_mix(ham_int, self.tracked_tips_50uL, self.fragmentation_buffer_positions, 
                mix_volume=200, mix_cycles=10,
                liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')

        # Multidispense reagents into MIDI plate stored on CPAC
        mix_position = [(self.MIDI_CPAC, idx) for idx in range(self.num_samples)]

        # Add Buffer EB (room temp)
        print("Adding Fragmentation Mix components...")
        multi_dispense(ham_int, tips=self.tracked_tips_300uL, source_positions=self.buffer_eb_positions, dispense_positions=mix_position,
                       volumes=[buffer_eb_vol]*self.num_samples, aspiration_height=0,
                       liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')

        # Add Fragmentation Buffer (room temp)
        multi_dispense(ham_int, tips=self.tracked_tips_300uL, source_positions=self.fragmentation_buffer_positions, dispense_positions=mix_position,
                       volumes=[frag_buffer_vol]*self.num_samples, aspiration_height=0,
                       liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
        
        # Add Fragmentation Enzyme (cold reagent from CPAC)
        multi_dispense(ham_int, self.tracked_tips_300uL, self.fragmentation_enzyme_positions, mix_position,
                       volumes=[frag_enzyme_vol]*self.num_samples, aspiration_height=0,
                       liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')

        # Mix the prepared fragmentation mix
        total_vol = buffer_eb_vol + frag_buffer_vol + frag_enzyme_vol
        pip_mix(ham_int, tips = self.tracked_tips_300uL,positions_to_mix=mix_position, mix_cycles=3, mix_volume=25,
                    liquid_height=0, liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')
        
        return mix_position

    def fragmentation_end_repair_atailing(self):
        """Step 1: GEX Fragmentation, End Repair & A-tailing."""
        if not self._initialized:
            print("Warning: Protocol not initialized. Call initialize() first.")
        
        print("Starting Fragmentation, End Repair & A-tailing...")
        
        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())

            # Prepare fragmentation mix
            self._prepare_fragmentation_mix(ham_int)
            
            # Transport MIDI plate from the TEC to MIDI_Pipette3
            transport_resource(ham_int, self.MIDI_CPAC, self.MIDI_Plate,
                               resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Transfer fragmentation mix to HSP_Plate2
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                         source_plate=self.MIDI_Plate, target_plate=self.HSP_Plate2, volume=40,
                         liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')

            # Move to thermal cycler
            odtc_open_door(ham_int, device_id=1)
            transport_resource(ham_int, self.HSP_Plate2, self.HSP_ODTC,
                              grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            
            # Add lid
            transport_resource(ham_int, self.Lid_Stack.fetch_next(), self.HSP_ODTC_Lid,
                              grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.LID, iswap=True)
            
            odtc_close_door(ham_int, device_id=1)
            
            # Run thermal cycling protocol: 32°C 5min, 65°C 30min
            odtc_execute_protocol(ham_int, device_id=1, method_name='odtc_protocols/10x_fragmentation_protocol.xml', simulating=self.device_simulation)

            # ODTC protocol is expected to run through to the beginning of the next step, we
            # multi-dispense reagents in the next step so we are ready to stamp them out as soon as
            # the ODTC protocol is done.
            
            print("Fragmentation, End Repair & A-tailing completed.")

    def post_fragmentation_spriselect(self):
        """Step 2: Post Fragmentation Double Sided SPRIselect."""
        print("Starting Post Fragmentation SPRIselect cleanup...")

        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())

            # Move HSP plate from stack to HSP_Plate
            transport_resource(ham_int, self.HSP_Stack.fetch_next(), self.HSP_Plate,
                               resource_type=GrippedResource.PCR, core_gripper=True)

            # Move MIDI plate from MIDI_Plate to HHS5_MIDI
            transport_resource(ham_int, self.MIDI_Plate, self.HHS5_MIDI,
                               resource_type=GrippedResource.MIDI, core_gripper=True)

            # Get MIDI plate from stack and move it to MIDI_Plate
            transport_resource(ham_int, self.MIDI_Stack.fetch_next(), self.MIDI_Plate,
                               resource_type=GrippedResource.MIDI, core_gripper=True)

            
            # Multidispense EB buffer to HSP_Plate3
            hsp_positions = [(self.MIDI_Plate, i) for i in range(self.num_samples)]
            multi_dispense(ham_int, self.tracked_tips_300uL, self.buffer_eb_positions, hsp_positions,
                           volumes=[40] * self.num_samples,
                           liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')

            
            # Move MIDI_Plate containing SPRI beads to HHS5_MIDI
            transport_resource(ham_int, self.MIDI_Plate, self.HHS5_MIDI,
                               resource_type=GrippedResource.MIDI, core_gripper=True)

            # Get MIDI plate from stack and move it to MIDI_Plate
            transport_resource(ham_int, self.MIDI_Stack.fetch_next(), self.MIDI_Plate,
                               resource_type=GrippedResource.MIDI, core_gripper=True)

            # Multidispense SPRI beads to MIDI_Plate (0.7x)
            midi_positions = [(self.MIDI_Plate, i) for i in range(self.num_samples)]
            multi_dispense(ham_int, self.tracked_tips_300uL, self.spriselect_positions, midi_positions,
                           volumes=[21] * self.num_samples,
                           liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')

            # Transport MIDI_Plate with beads to CPAC while we wait for ODTC protocol to finish
            transport_resource(ham_int, self.MIDI_Plate, self.MIDI_CPAC,
                               resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Wait for ODTC protocol to finish
            odtc_wait_for_idle(ham_int, device_id=1, simulating=self.device_simulation) # Run dehumidify protocol?

            # Remove from thermal cycler
            odtc_open_door(ham_int, device_id=1)

            transport_resource(ham_int, self.HSP_ODTC_Lid, self.Lid_Stack.put_back(),
                               grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.LID, iswap=True)
            transport_resource(ham_int, self.HSP_ODTC, self.HSP_Plate2,
                               grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            
            odtc_close_door(ham_int, device_id=1)

            # Transport MIDI plate from CPAC_03 to MIDI_Plate
            transport_resource(ham_int, self.MIDI_CPAC, self.MIDI_Plate,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Mix MIDI_Plate
            mix_plate(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                      self.MIDI_Plate, mixing_volume=75, mix_cycles=15,
                      liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
            
            # Transfer 30 uL pre-plated beads from MIDI_Plate to HSP_Plate2
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                        self.MIDI_Plate, self.HSP_Plate2, volume=30, dispense_mix_cycles=15,
                        dispense_mix_volume=75, liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')

            # Incubate 5 minutes
            incubate_timer = start_timer(300)
            incubate_timer.wait(skip=self.device_simulation)

            # Move HSP_Plate2 to HSP_Magnet
            transport_resource(ham_int, self.HSP_Plate2, self.HSP_Magnet,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Wait for beads to settle
            settle_timer = start_timer(180)
            settle_timer.wait(skip=self.device_simulation)

            # Remove supernatant (70μl) and transfer to HSP_Plate
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                        self.HSP_Magnet, self.HSP_Plate, volume=70, liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')

            # Mix beads (15x at 75μl)
            mix_plate(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                      self.MIDI_Plate, mixing_volume=75, mix_cycles=15,
                      liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')

            # Transfer second round of beads (10uL)
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                        self.MIDI_Plate, self.HSP_Plate2, volume=10, dispense_mix_cycles=15, dispense_mix_volume=30,
                        liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')

            # Incubate 5 minutes
            incubate_timer = start_timer(300)
            incubate_timer.wait(skip=self.device_simulation)
            
            # Move HSP_Plate2 to HSP_Magnet
            transport_resource(ham_int, self.HSP_Plate2, self.HSP_Magnet,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Wait for beads to settle
            settle_timer = start_timer(120)

            # Transport MIDI_Plate to CPAC_03
            transport_resource(ham_int, self.MIDI_Plate, self.MIDI_CPAC,
                               resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Transport MIDI plate from HHS5_MIDI to MIDI_Plate
            transport_resource(ham_int, self.HHS5_MIDI, self.MIDI_Plate,
                               resource_type=GrippedResource.MIDI, core_gripper=True)
            # Incubate 3 minutes
            settle_timer.wait(skip=self.device_simulation)
            
            # Remove supernatant (150μl) and transfer to HSP_Plate
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                        self.HSP_Magnet, self.Liquid_Waste, volume=150, liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')

            # Ethanol washes (2x with 125μl)
            for wash in range(2):
                transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                            self.EthanolReservoir, self.HSP_Magnet, volume=125,
                            liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
                
                # Wait 30 seconds
                wait_timer = start_timer(30)
                wait_timer.wait(skip=self.device_simulation)

                # Remove ethanol
                double_aspirate_supernatant_96(ham_int, tips=self.tracked_tips_300uL, tip_support=self.tip_support, num_samples=self.num_samples,
                                              source_plate=self.HSP_Magnet, destination_plate=self.Liquid_Waste, first_volume=100, second_volume=100,
                                              liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            # Air dry briefly (don't over-dry)
            dry_timer = start_timer(60)
            dry_timer.wait(skip=self.device_simulation)

            # Move off magnet
            transport_resource(ham_int, self.HSP_Magnet, self.HSP_Plate2,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            

            # Add EB buffer (32μl)
            sample_positions = [(self.HSP_Plate2, i) for i in range(self.num_samples)]
            pip_transfer(ham_int, self.tracked_tips_300uL, self.buffer_eb_positions, sample_positions, volumes=[32]*self.num_samples,
                        liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            # Transport HSP_Plate2 to HHS1_HSP
            transport_resource(ham_int, self.HSP_Plate2, self.HHS1_HSP,
                               resource_type=GrippedResource.PCR, core_gripper=True)

            # Start incubation timer
            incubate_timer = start_timer(120)

            # Get HSP plate from stack and move it to HSP_Pipette
            transport_resource(ham_int, self.HSP_Stack.fetch_next(), self.HSP_Plate,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Wait for incubation to complete
            incubate_timer.wait(skip=self.device_simulation)

            # Transport HSP_Plate from HHS1 to HSP_Magnet
            transport_resource(ham_int, self.HHS1_HSP, self.HSP_Magnet,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Wait for beads to settle
            settle_timer = start_timer(180)
            settle_timer.wait(skip=self.device_simulation)

            # Transfer supernatant (32μl) to HSP_Plate
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                        self.HSP_Magnet, self.HSP_Plate, volume=32, liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')
            print("Post Fragmentation SPRIselect cleanup completed.")

    def _prepare_ligation_mix(self, ham_int):
        """Prepare adapter ligation mix with 10% excess."""
        total_reactions = int(self.num_samples * self.excess_factor)
        
        ligation_mix_vol = 40 * total_reactions
        dna_ligase_vol = 10 * total_reactions
        
        # Prepare mix in position 4 of CAR_VIALS_SMALL
        mix_position = [(self.CAR_VIALS_SMALL, 4)]
        
        # Add Ligation Mix (room temp)
        pip_transfer(ham_int, self.tracked_tips_1000uL, self.ligation_mix_positions, mix_position,
                     volumes=[ligation_mix_vol], liquid_class='HighVolumeFilter_Water_DispenseSurface_Empty')
        
        # Add DNA Ligase (cold reagent from CPAC)
        pip_transfer(ham_int, self.tracked_tips_1000uL, self.dna_ligase_positions, mix_position,
                     volumes=[dna_ligase_vol], liquid_class='HighVolumeFilter_Water_DispenseSurface_Empty')

        # Mix
        total_vol = ligation_mix_vol + dna_ligase_vol
        pip_mix(ham_int, self.tracked_tips_1000uL, mix_position,
                mix_volume=min(total_vol, 1000), mix_cycles=15,
                liquid_class='HighVolumeFilter_Water_DispenseSurface_Empty')

        return mix_position

    def adapter_ligation(self):
        """Step 3: GEX Adapter Ligation."""
        print("Starting Adapter Ligation...")

        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())
            
            # Move HSP_Plate to HHS1_HSP
            transport_resource(ham_int, self.HSP_Plate, self.HHS1_HSP,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Get MIDI plate from stack and move it to MIDI_CPAC
            transport_resource(ham_int, self.MIDI_Stack.fetch_next(), self.MIDI_CPAC,
                               resource_type=GrippedResource.MIDI, core_gripper=True)


            # Transport MIDI plate from CPAC_02 to MIDI_Pipette3 because otherwise it will be in the way of HSP ODTC
            transport_resource(ham_int, self.MIDI_CPAC, self.MIDI_Plate,
                               resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Prepare ligation mix
            ligation_mix_pos = self._prepare_ligation_mix(ham_int)

            # Channels transfer from ligation mix to HSP_Plate (50μl)
            HSP_Plate_positions = [(self.HSP_Plate, i) for i in range(self.num_samples)]
            pip_transfer(ham_int, self.tracked_tips_50uL, ligation_mix_pos, HSP_Plate_positions, volumes=[50]*self.num_samples, 
                         liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')

            # Mix (15x at 90μl)
            mix_plate(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                      self.HSP_Plate, mixing_volume=90, mix_cycles=15,
                      liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')


            # Move to thermal cycler
            odtc_open_door(ham_int, device_id=1)
            transport_resource(ham_int, self.HSP_Plate2, self.HSP_ODTC,
                              grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            
            # Add lid
            transport_resource(ham_int, self.Lid_Stack.fetch_next(), self.HSP_ODTC_Lid,
                              grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.LID, iswap=True)
            
            odtc_close_door(ham_int, device_id=1)
            
            # Run ligation protocol (20°C for 15 min)
            odtc_execute_protocol(ham_int, device_id=1, method_name='odtc_protocols/10x_adapter_ligation_protocol.xml', simulating=self.device_simulation)

            # Remove from thermal cycler
            odtc_open_door(ham_int, device_id=1)
            transport_resource(ham_int, self.HSP_ODTC_Lid, self.Lid_Stack.put_back(),
                              grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.LID, iswap=True)
            transport_resource(ham_int, self.HSP_ODTC, self.HSP_Plate2,
                              grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            odtc_close_door(ham_int, device_id=1)

            # Transport MIDI plate from MIDI_Pipette3 to MIDI_CPAC
            transport_resource(ham_int, self.MIDI_Plate, self.MIDI_CPAC,
                               resource_type=GrippedResource.PCR, core_gripper=True)

            print("Adapter Ligation completed.")

    def post_ligation_cleanup(self):
        """Step 4: Post Ligation Cleanup - SPRIselect."""
        print("Starting Post Ligation SPRIselect cleanup...")

        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())

            # Transfer samples (100μl) to MIDI plate
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                        self.HSP_Plate, self.MIDI_Plate, volume=100,
                        liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
            
            # multidispense EB buffer to MIDI plate 3
            midi_positions = [(self.MIDI_Plate, i) for i in range(self.num_samples)]
            multi_dispense(ham_int, self.tracked_tips_300uL, self.buffer_eb_positions, midi_positions,
                           volumes=[30.5] * self.num_samples, liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
            
            # move MIDI plate 3 to HHS5_MIDI
            transport_resource(ham_int, self.MIDI_Plate, self.HHS5_MIDI,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Get MIDI plate from stack and move it to MIDI_Plate
            self.MIDI_Stack.reset_all()
            transport_resource(ham_int, self.MIDI_Stack.fetch_next(), self.MIDI_Plate,
                               resource_type=GrippedResource.MIDI, core_gripper=True)

            # Add SPRIselect (0.8X = 80μl)
            midi_positions = [(self.MIDI_Plate, i) for i in range(self.num_samples)]
            multi_dispense(ham_int, self.tracked_tips_300uL, self.spriselect_positions, midi_positions,
                           volumes=[80] * self.num_samples,
                           liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
            
            # Mix (15x at 150μl)
            mix_plate(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                      self.MIDI_Plate, mixing_volume=150, mix_cycles=15,
                      liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
                        

            # Transport MIDI plate from CPAC to MIDI_Plate
            transport_resource(ham_int, self.MIDI_CPAC, self.MIDI_Plate,
                               resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Mix plate with beads (15x at 150μl)
            mix_plate(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                      self.MIDI_Plate, mixing_volume=150, mix_cycles=15,
                      liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
            
            # Transfer beads to HSP plate
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                        self.MIDI_Plate, self.HSP_Plate2, volume=150, liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
            
            # Transport HSP_Plate2 to HSP_Magnet
            transport_resource(ham_int, self.HSP_Plate2, self.HSP_Magnet,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Incubate 4 minutes
            incubate_timer = start_timer(240)

            # Transport MIDI Plate 3 to MIDI CPAC
            transport_resource(ham_int, self.MIDI_Plate, self.MIDI_CPAC,
                               resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Move MIDI plate from HHS5 to MIDI Plate 3
            transport_resource(ham_int, self.HHS5_MIDI, self.MIDI_Plate,
                               resource_type=GrippedResource.MIDI, core_gripper=True)
            
            incubate_timer.wait(skip=self.device_simulation)
            
            # Double aspirate supernatant from HSP_Magnet
            double_aspirate_supernatant_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                                          self.HSP_Magnet, self.MIDI_Waste, first_volume=100, second_volume=100,
                                          liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            
            # Ethanol washes (2x with 180μl)
            magnet_positions = [(self.HSP_Magnet, i) for i in range(self.num_samples)]
            for wash in range(2):
                transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                            self.EthanolReservoir, self.HSP_Magnet, volume=180,
                            liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
                
                wait_timer = start_timer(30)
                wait_timer.wait(skip=self.device_simulation)
                
                transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                            self.HSP_Magnet, self.MIDI_Waste, volume=180,
                            liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                            aspiration_height=0.3)
            
            # Drying timer
            settle_timer = start_timer(120)
            settle_timer.wait(skip=self.device_simulation)

            # Move off magnet
            transport_resource(ham_int, self.HSP_Magnet, self.HSP_Plate2,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Add Buffer EB (30.5μl)
            offmagnet_positions = [(self.HSP_Plate2, i) for i in range(self.num_samples)]
            pip_transfer(ham_int, self.tracked_tips_300uL, self.buffer_eb_positions, offmagnet_positions,
                           volumes=[30.5] * self.num_samples, mix_cycles=15, vol_mix_dispense=25,
                           liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')
            
            # Incubate 2 minutes
            incubate_timer = start_timer(120)
            incubate_timer.wait(skip=self.device_simulation)

            # Move to magnet
            transport_resource(ham_int, self.HSP_Plate2, self.HSP_Magnet,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Wait for beads to settle
            settle_timer = start_timer(240)

            # Transport MIDI plate from MIDI_Plate to MIDI_HHS5
            transport_resource(ham_int, self.MIDI_Plate, self.HHS5_MIDI,
                               resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Get HSP plate from stack and move it to HSP_Pipette2
            transport_resource(ham_int, self.HSP_Stack.fetch_next(), self.HSP_Plate2,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            

            settle_timer.wait(skip=self.device_simulation)

            # Remove supernatant from HSP_Magnet to HSP_Plate2
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                        self.HSP_Magnet, self.HSP_Plate2, volume=30, liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')
            
            print("Post Ligation SPRIselect cleanup completed.")

    def sample_index_pcr(self):
        """Step 5: Sample Index PCR."""
        print("Starting Sample Index PCR...")

        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())

            # Move HSP_Plate2 to HHS1_HSP
            transport_resource(ham_int, self.HSP_Plate2, self.HHS1_HSP,
                               resource_type=GrippedResource.PCR, core_gripper=True)

            # Get MIDI plate from stack and move it to MIDI_CPAC
            transport_resource(ham_int, self.MIDI_Stack.fetch_next(), self.MIDI_CPAC,
                               resource_type=GrippedResource.MIDI, core_gripper=True)

            # Multidispense Library Amp Mix (50μl) to intermediate HSP plate
            intermediate_positions = [(self.HSP_Plate2, i) for i in range(self.num_samples)]
            multi_dispense(ham_int, self.tracked_tips_300uL, self.library_amp_mix_positions,
                            intermediate_positions, volumes=[self.amp_mix_volume] * self.num_samples,
                            liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')

            # Move MIDI plate from CPAC to MIDI_Pipette3
            transport_resource(ham_int, self.MIDI_CPAC, self.MIDI_Plate,
                               resource_type=GrippedResource.MIDI, core_gripper=True)

            # Pip transfer from index plate on HHS1 to sample plate
            pip_transfer(ham_int, self.tracked_tips_50uL, self.index_positions, intermediate_positions,
                         volumes=[self.index_volume] * self.num_samples,
                         liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')

            
            # Move to thermal cycler
            odtc_open_door(ham_int, device_id=1)
            transport_resource(ham_int, self.HSP_Plate, self.HSP_ODTC,
                              grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            
            # Add lid
            transport_resource(ham_int, self.Lid_Stack.fetch_next(), self.HSP_ODTC_Lid,
                              grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.LID, iswap=True)
            
            odtc_close_door(ham_int, device_id=1)
            
            # Run PCR protocol (cycles based on cDNA input)
            odtc_execute_protocol(ham_int, device_id=1, method_name='odtc_protocols/10x_sample_index_pcr.xml', simulating=self.device_simulation)

            # Remove from thermal cycler
            odtc_open_door(ham_int, device_id=1)
            transport_resource(ham_int, self.HSP_ODTC_Lid, self.Lid_Stack.put_back(),
                              grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.LID, iswap=True)
            transport_resource(ham_int, self.HSP_ODTC, self.HSP_Plate2,
                              grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
            odtc_close_door(ham_int, device_id=1)

            print("Sample Index PCR completed.")

    def final_size_selection(self):
        """Step 6: Post Sample Index PCR Double Sided Size Selection."""
        print("Starting Final Size Selection...")

        with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
            ham_int.initialize()
            normal_logging(ham_int, os.getcwd())

            # Move HSP_Plate to HHS1_HSP
            transport_resource(ham_int, self.HSP_Plate2, self.HHS1_HSP,
                               resource_type=GrippedResource.PCR, core_gripper=True)

            # Get HSP plate from stack and move it to HSP_Plate
            transport_resource(ham_int, self.HSP_Stack.fetch_next(), self.HSP_Plate,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Multi-dispense EB buffer to MIDI plate 3
            midi_positions = [(self.MIDI_Plate, i) for i in range(self.num_samples)]
            multi_dispense(ham_int, self.tracked_tips_300uL, self.buffer_eb_positions, midi_positions,
                           volumes=[35.5] * self.num_samples, liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
            
            # Move MIDI plate 3 to HHS5_MIDI
            transport_resource(ham_int, self.MIDI_Plate, self.HHS5_MIDI,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            # Get MIDI plate from stack and move it to MIDI_Plate
            transport_resource(ham_int, self.MIDI_Stack.fetch_next(), self.MIDI_Plate,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Multi-dispense SPRI beads to MIDI_Plate (0.6x)
            midi_positions = [(self.MIDI_Plate, i) for i in range(self.num_samples)]
            multi_dispense(ham_int, self.tracked_tips_300uL, self.spriselect_positions, midi_positions,
                           volumes=[60] * self.num_samples, liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
            

            # Transport MIDI_Plate to CPAC_03
            transport_resource(ham_int, self.MIDI_Plate, self.MIDI_CPAC,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Wait for ODTC protocol to complete
            odtc_ready = False
            while not odtc_ready:
                status = odtc_get_status(ham_int, device_id=1, simulating=self.device_simulation)
                if status.state == 'idle':
                    odtc_ready = True
                else:
                    time.sleep(5)

            # Move MIDI plate from CPAC_02 to MIDI_Pipette3
            transport_resource(ham_int, self.MIDI_CPAC, self.MIDI_Plate,
                               resource_type=GrippedResource.MIDI, core_gripper=True)

            # Mix SPRI beads in MIDI plate (15x at 150μl)
            mix_plate(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                      self.MIDI_Plate, mixing_volume=150, mix_cycles=15, liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            # Transfer beads to HSP plate
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                        self.MIDI_Plate, self.HSP_Plate2, volume=150,
                        liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            # Incubate 5 minutes
            incubate_timer = start_timer(300)
            incubate_timer.wait(skip=self.device_simulation)

            # Move HSP plate to magnet
            transport_resource(ham_int, self.HSP_Plate2, self.HSP_Magnet,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Wait for beads to settle
            settle_timer = start_timer(180)
            settle_timer.wait(skip=self.device_simulation)

            # Aspirate supernatant (150μl) and transfer to HSP_Plate
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                        self.HSP_Magnet, self.HSP_Plate, volume=150, liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')


            # Transport HSP_Magnet to HSP_Waste
            transport_resource(ham_int, self.HSP_Magnet, self.HSP_Waste,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            
            
            # Add second SPRI beads (0.2x = 20μl)
            transfer_96(ham_int, self.tracked_tips_300uL,self.tip_support, self.num_samples,
                        self.MIDI_Plate, self.HSP_Plate,
                        volume=20, dispense_mix_cycles=10, dispense_mix_volume=15,
                        liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            # Incubate beads 2 minutes
            incubate_timer = start_timer(120)
            incubate_timer.wait(skip=self.device_simulation)

            # Move HSP_Plate to HSP_Magnet
            transport_resource(ham_int, self.HSP_Plate, self.HSP_Magnet,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Wait for beads to settle
            settle_timer = start_timer(180)
            settle_timer.wait(skip=self.device_simulation)

            # Transfer supernatant (150μl) to waste
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                        self.HSP_Magnet, self.HSP_Waste, volume=150,
                        liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')
            
            # Ethanol washes (2x with 200μl)
            magnet_positions = [(self.HSP_Magnet, i) for i in range(self.num_samples)]
            for wash in range(2):
                transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                            self.EthanolReservoir, self.HSP_Magnet, volume=200,
                            liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')
                
                wait_timer = start_timer(30)
                wait_timer.wait(skip=self.device_simulation)
                
                transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                            self.HSP_Magnet, self.MIDI_Waste, volume=200,
                            liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                            aspiration_height=0.3)
            
            # Drying timer
            drying_timer = start_timer(120)
            drying_timer.wait(skip=self.device_simulation)

            # Move off magnet
            transport_resource(ham_int, self.HSP_Magnet, self.HSP_Plate2,
                               resource_type=GrippedResource.PCR, core_gripper=True)

            # Move MIDI plate from MIDI_Pipette3 to MIDI_CPAC
            transport_resource(ham_int, self.MIDI_Plate, self.MIDI_CPAC,
                               resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Move MIDI plate from HHS5 to MIDI_Plate
            transport_resource(ham_int, self.HHS5_MIDI, self.MIDI_Plate,
                               resource_type=GrippedResource.MIDI, core_gripper=True)
            
            settle_timer.wait(skip=self.device_simulation)

            # Add Buffer EB (35.5μl)
            transfer_96(ham_int, self.tracked_tips_300uL,self.tip_support, self.num_samples,
                        self.MIDI_Plate, self.HSP_Plate2,
                        volume=35.5, dispense_mix_cycles=15, dispense_mix_volume=25,
                        liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            # Incubate 2 minutes
            incubate_timer = start_timer(120)
            incubate_timer.wait(skip=self.device_simulation)

            # Move MIDI plate 3 back to HHS5_MIDI
            transport_resource(ham_int, self.MIDI_Plate, self.HHS5_MIDI,
                               resource_type=GrippedResource.MIDI, core_gripper=True)
            
            # Move to magnet
            transport_resource(ham_int, self.HSP_Plate2, self.HSP_Magnet,
                               resource_type=GrippedResource.PCR, core_gripper=True)
            
            # Wait for beads to settle
            settle_timer = start_timer(240)

            # Get HSP plate from stack and move it to HSP_Plate2
            transport_resource(ham_int, self.HSP_Stack.fetch_next(), self.HSP_Plate2,
                               resource_type=GrippedResource.PCR, core_gripper=True)

            settle_timer.wait(skip=self.device_simulation)


            # Remove supernatant (165μl)
            transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                        self.HSP_Magnet, self.HSP_Plate2, volume=165,
                        liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            print("Final libraries ready for QC and sequencing!")
            print("Final Size Selection completed.")

    def run_complete_protocol(self):
        """Run the complete 10X GEX library construction protocol."""
        print(f"Starting 10X GEX Library Construction for {self.num_samples} samples")
        print(f"Sample volume: {self.sample_volume} µL")
        print(f"Using {self.pcr_cycles} PCR cycles")
        print(f"Simulation mode: {self.simulation}")
        
        try:
            # Run protocol steps
            self.initialize()
            self.fragmentation_end_repair_atailing()
            self.post_fragmentation_spriselect()
            self.adapter_ligation()
            self.post_ligation_cleanup()
            self.sample_index_pcr()
            self.final_size_selection()
            
            print("10X GEX Library Construction Protocol Completed Successfully!")
            print("Libraries are ready for quality control and sequencing.")
            
        except Exception as e:
            print(f"Protocol failed with error: {e}")
            raise


# Example usage
if __name__ == "__main__":
    # Initialize protocol with 8 samples and 12 PCR cycles
    # PCR cycles should be optimized based on cDNA input:
    # 0.25-50 ng: 14-16 cycles
    # 50-250 ng: 12-14 cycles  
    # 250-600 ng: 10-12 cycles
    # 600-1,100 ng: 8-10 cycles
    # 1,100-1,500 ng: 6-8 cycles
    # >1,500 ng: 5 cycles
    
    protocol = TenXGEXLibraryPrepProtocol(
        num_samples=8,
        sample_volume=10,
        pcr_cycles=12,
        simulation=True,
        device_simulation=True
    )

    protocol.run_protocol()
    generate_tadm_report()