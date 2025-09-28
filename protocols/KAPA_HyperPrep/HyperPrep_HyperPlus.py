# KAPA HyperPlus Library Preparation Protocol

from pyhamilton import (HamiltonInterface, LayoutManager, Plate96, Tip96, 
                       TrackedTips, StackedResources, normal_logging, layout_item, start_timer)

from pyhamilton.devices import (initialize_cpac, set_temperature_target_cpac, start_temperature_control_cpac,
                                 get_temperature_cpac, odtc_initialize, odtc_connect, odtc_open_door, odtc_close_door, odtc_execute_protocol,
                                 odtc_wait_for_idle)

from pyhamilton.consumables import (tracked_volume_aspirate, ReagentTrackedReservoir60mL, ReagentTrackedPlate96, 
                                    ReagentTrackedFalconCarrier24, generate_reagent_summary, ReagentTrackedBulkPlate, ReagentTrackedEppiCarrier32)
from pyhamilton.ngs import Protocol, LoadingVis, generate_tadm_report

from pyhamilton.resources import TrackedTips, StackedResources, Plate96, Tip96, Waste96, LayoutManager,TipSupportTracker, Lid, layout_item

from pyhamilton.pipetting import (pip_transfer, multi_dispense, double_aspirate_supernatant_96, ethanol_wash, transfer_96, pip_mix, mix_plate,
                                  shear_plate_96, transfer_96)
from pyhamilton.transport import (transport_resource, GripDirection, GrippedResource, GripperParams)

import time
import os


class KAPAHyperPlusProtocol(Protocol):
   """
   KAPA HyperPlus Library Preparation protocol class for Hamilton liquid handler.
   
   This class encapsulates the complete workflow for KAPA HyperPlus library preparation
   including enzymatic fragmentation, end repair & A-tailing, adapter ligation, 
   post-ligation cleanup, library amplification, and final cleanup.
   """
   
   def __init__(self, num_samples=96, dna_input_ng=100, fragmentation_time_min=15, pcr_cycles=8, simulation=False, device_simulation=False):
       """
       Initialize the protocol with sample parameters.
       
       Args:
           num_samples (int): Number of samples to process (default: 96)
           dna_input_ng (int): DNA input amount in ng (1-1000 ng range) (default: 100)
           fragmentation_time_min (int): Fragmentation time in minutes for target insert size (default: 15)
           pcr_cycles (int): Number of PCR cycles based on input amount (default: 8)
           simulation (bool): Whether to run in simulation mode (default: False)
       """

       self.available_steps = [
           ("Initialize System", "initialize"),
           ("Enzymatic Fragmentation", "enzymatic_fragmentation"),
           ("End Repair & A-Tailing", "end_repair_a_tailing"),
           ("Adapter Ligation", "adapter_ligation"),
           ("Post-Ligation Cleanup", "post_ligation_cleanup"),
           ("Library Amplification", "library_amplification"),
           ("Final Cleanup & Size Selection", "final_cleanup_size_selection")
       ]

       super().__init__()
       self.num_samples = num_samples
       self.dna_input_ng = dna_input_ng
       self.fragmentation_time_min = fragmentation_time_min
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
       self.lmgr = LayoutManager('KAPA_HyperPlus_Deck.lay')

       # Deck Resources - using only allowed positions
       self.MIDI_OnMagnet = layout_item(self.lmgr, Plate96, 'MIDI_OnMagnet')
       self.MIDI_OffMagnet = layout_item(self.lmgr, Plate96, 'MIDI_Pipette')
       self.HSP_Plate = layout_item(self.lmgr, Plate96, 'HSP_Pipette')
       self.HSP_Plate_2 = layout_item(self.lmgr, Plate96, 'HSP_Pipette2')

       self.HSP_ODTC = layout_item(self.lmgr, Plate96, 'HSP_ODTC')
       self.ODTC_Lid = layout_item(self.lmgr, Lid, 'Ham_ComfortLid_ODTC')

       # Reagent containers with proper naming
       self.CAR_VIALS_SMALL = layout_item(self.lmgr, ReagentTrackedEppiCarrier32, 'SMP_CAR_24_15x75_A00_0001') # Use this for small volume reagents
       self.CPAC_Reagents = layout_item(self.lmgr, ReagentTrackedPlate96, 'CPAC_HSP_0001') # Use this plate for cold reagents
       self.EthanolReservoir = layout_item(self.lmgr, ReagentTrackedBulkPlate, 'RGT_Ethanol') # Use this plate for ethanol
       self.BeadMix_Container = layout_item(self.lmgr, ReagentTrackedReservoir60mL, 'rgt_cont_60ml_BC_A00_0002')

       # Reagent mapping and positions
       # Pre-mixed reagents (enzymes and buffers combined manually before protocol)
       self.fragmentation_mastermix_positions = self.CAR_VIALS_SMALL.assign_reagent_map('FragmentationMasterMix', [0])
       self.end_repair_mastermix_positions = self.CAR_VIALS_SMALL.assign_reagent_map('EndRepairMasterMix', [1])
       self.ligation_mastermix_positions = self.CAR_VIALS_SMALL.assign_reagent_map('LigationMasterMix', [2])
       self.kapa_adapters_positions = self.CAR_VIALS_SMALL.assign_reagent_map('KAPA_Adapters', [3])
       self.nuclease_free_water_positions = self.CAR_VIALS_SMALL.assign_reagent_map('NucleaseFreeWater', [4])
       
       # Assign index primers to an untracked HSP plate
       self.index_primer_positions = [(self.HSP_Plate_2, i) for i in range(self.num_samples)]


       # Cold-sensitive reagents in CPAC (4°C)
       self.kapa_hifi_mix_positions = self.CPAC_Reagents.assign_reagent_map('KAPA_HiFi_Mix', [0])

       # Bulk reagents with proper trough naming
       self.ethanol_plate = self.EthanolReservoir.assign_reagent_map('Ethanol80', range(96))
       self.kapa_pure_beads_positions = self.BeadMix_Container.assign_reagent_map('KAPA_Pure_Beads', range(8))


       # Stacked resources
       self.Lid_Stack = StackedResources.from_prefix(
           tracker_id="Ham_ComfortLid_Stack", 
           prefix="Ham_ComfortLid_Stack",
           lmgr=self.lmgr,
           resource_type=Lid,
           count=4)

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

       # Tip support with correct naming
       self.tip_support_resource = layout_item(self.lmgr, Tip96, 'TipSupport_0001')
       self.tip_support = TipSupportTracker(self.tip_support_resource)

       # Batch together tracked objects for resource consumption estimates
       self.tracked_reagent_vessels = [
           self.CAR_VIALS_SMALL, self.CPAC_Reagents, self.EthanolReservoir, 
           self.BeadMix_Container
       ]
       self.tracked_tips = [self.tracked_tips_50uL, self.tracked_tips_300uL]
       self.stacked_resources = [self.Lid_Stack]

   def _setup_volumes(self):
       """Setup volume parameters based on input and protocol requirements."""
       # DNA input volume (assuming 1ng/µL concentration for calculation)
       self.dna_input_volume = min(self.dna_input_ng, 50)  # Max 50µL input
       
       # Reagent volumes per sample (pre-mixed master mixes)
       self.fragmentation_mastermix_volume = 10  # µL (buffer + enzyme pre-mixed)
       self.end_repair_mastermix_volume = 10     # µL (buffer + enzyme pre-mixed)
       self.ligation_mastermix_volume = 20       # µL (buffer + enzyme pre-mixed)
       self.adapter_volume = 2.5                 # µL
       self.pcr_mix_volume = 25                  # µL
       self.index_primer_volume = 5              # µL
       
       # Bead volumes for cleanups
       self.post_ligation_bead_volume = 50       # 1.0X ratio
       self.final_cleanup_bead_volume = 45       # 0.9X ratio for size selection
       self.ethanol_wash_volume = 200            # µL per wash
       
       # Elution volumes
       self.post_ligation_elution_volume = 25    # µL
       self.final_elution_volume = 22            # µL
       
       # Calculated volumes with excess
       self.excess_factor = 1.1

   def initialize(self):
       """Initialize Hamilton system and CPAC for cold reagents."""
       with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
           ham_int.initialize()
           normal_logging(ham_int, os.getcwd())
           
           # Initialize CPAC for cold reagent storage
           initialize_cpac(ham_int, controller_id=1, simulating=self.device_simulation)
           set_temperature_target_cpac(ham_int, controller_id=1, device_id=1, target_temp=4)
           start_temperature_control_cpac(ham_int, controller_id=1, device_id=1)
           
           # Verify CPAC temperature
           temp = get_temperature_cpac(ham_int, controller_id=1, device_id=1)
           print(f"CPAC temperature stabilized at: {temp}°C")

           device_id = odtc_connect(ham_int, simulation_mode=self.device_simulation, local_ip='1.2.3.4', device_ip='5.6.7.8', device_port='COM4')  # Replace with actual IPs
           odtc_initialize(ham_int, device_id=device_id)
           
           self._initialized = True
           print("KAPA HyperPlus Protocol initialization completed successfully.")

   def enzymatic_fragmentation(self):
       """Step 1: Enzymatic Fragmentation."""
       if not self._initialized:
           print("Warning: Protocol not initialized. Call initialize() first.")
       
       print("Starting Enzymatic Fragmentation step...")
       
       with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, persistent=self.persistent) as ham_int:
           ham_int.initialize()
           normal_logging(ham_int, os.getcwd())

           # Transfer DNA samples from HSP_Plate (samples start here)
           sample_positions = [(self.HSP_Plate, i) for i in range(self.num_samples)]
           
           # Add Fragmentation Master Mix (pre-mixed buffer + enzyme)
           pip_transfer(ham_int, self.tracked_tips_50uL, self.fragmentation_mastermix_positions, sample_positions,
                        volumes=[self.fragmentation_mastermix_volume] * self.num_samples,
                        liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')
           
           # Mix samples
           mix_plate(ham_int, self.tracked_tips_50uL, self.tip_support, self.num_samples,
                     self.HSP_Plate, mixing_volume=15, mix_cycles=10,
                     liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')
           
           # Move to thermal cycler for fragmentation
           odtc_open_door(ham_int, device_id=1)
           transport_resource(ham_int, self.HSP_Plate, self.HSP_ODTC,
                             resource_type=GrippedResource.PCR, core_gripper=True)
           
           # Add lid
           transport_resource(ham_int, self.Lid_Stack.fetch_next(), self.ODTC_Lid,
                             resource_type=GrippedResource.LID, core_gripper=True)
           
           odtc_close_door(ham_int, device_id=1)
           
           # Run fragmentation protocol (37°C for specified time)
           odtc_execute_protocol(ham_int, device_id=1, method_name='Enzymatic_Fragmentation', simulating=self.device_simulation)
           
           # Remove from thermal cycler
           odtc_open_door(ham_int, device_id=1)
           transport_resource(ham_int, self.ODTC_Lid, self.Lid_Stack.put_back(),
                             resource_type=GrippedResource.LID, core_gripper=True)
           transport_resource(ham_int, self.HSP_ODTC, self.HSP_Plate,
                             resource_type=GrippedResource.PCR, core_gripper=True)
           odtc_close_door(ham_int, device_id=1)

           print("Enzymatic Fragmentation completed.")

   def end_repair_a_tailing(self):
       """Step 2: End Repair & A-Tailing (combined step)."""
       print("Starting End Repair & A-Tailing step...")

       with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, persistent=self.persistent) as ham_int:
           ham_int.initialize()
           normal_logging(ham_int, os.getcwd())

           # Transfer fragmented DNA to HSP_Plate_2
           current_volume = self.dna_input_volume + self.fragmentation_mastermix_volume
           
           transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                       self.HSP_Plate, self.HSP_Plate_2, volume=current_volume,
                       liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')
           
           # Add End Repair Master Mix (pre-mixed buffer + enzyme)
           end_repair_positions = [(self.HSP_Plate_2, i) for i in range(self.num_samples)]
           pip_transfer(ham_int, self.tracked_tips_50uL, self.end_repair_mastermix_positions, end_repair_positions,
                        volumes=[self.end_repair_mastermix_volume] * self.num_samples,
                        liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')
           
           # Mix samples
           mix_plate(ham_int, self.tracked_tips_50uL, self.tip_support, self.num_samples,
                     self.HSP_Plate_2, mixing_volume=25, mix_cycles=10,
                     liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')
           
           # Move to thermal cycler for end repair & A-tailing
           odtc_open_door(ham_int, device_id=1)
           transport_resource(ham_int, self.HSP_Plate_2, self.HSP_ODTC,
                             resource_type=GrippedResource.PCR, core_gripper=True)
           
           # Add lid
           transport_resource(ham_int, self.Lid_Stack.fetch_next(), self.ODTC_Lid,
                             resource_type=GrippedResource.LID, core_gripper=True)
           
           odtc_close_door(ham_int, device_id=1)
           
           # Run end repair & A-tailing protocol (20°C 30min, 65°C 30min)
           odtc_execute_protocol(ham_int, device_id=1, method_name='End_Repair_A_Tailing', simulating=self.device_simulation)
           odtc_wait_for_idle(ham_int, device_id=1, check_interval=30, max_wait=7200, simulating=self.device_simulation)

           # Remove from thermal cycler
           odtc_open_door(ham_int, device_id=1)
           transport_resource(ham_int, self.ODTC_Lid, self.Lid_Stack.put_back(),
                             resource_type=GrippedResource.LID, core_gripper=True)
           transport_resource(ham_int, self.HSP_ODTC, self.HSP_Plate_2,
                             resource_type=GrippedResource.PCR, core_gripper=True)
           odtc_close_door(ham_int, device_id=1)

           print("End Repair & A-Tailing completed.")

   def adapter_ligation(self):
       """Step 3: Adapter Ligation."""
       print("Starting Adapter Ligation step...")

       with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
           ham_int.initialize()
           normal_logging(ham_int, os.getcwd())

           # Transfer end-repaired DNA back to HSP_Plate
           current_volume = (self.dna_input_volume + self.fragmentation_mastermix_volume + 
                           self.end_repair_mastermix_volume)
           
           transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                       self.HSP_Plate_2, self.HSP_Plate, volume=current_volume,
                       liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')
           
           # Add KAPA Adapters
           ligation_positions = [(self.HSP_Plate, i) for i in range(self.num_samples)]
           pip_transfer(ham_int, self.tracked_tips_50uL, self.kapa_adapters_positions, ligation_positions,
                        volumes=[self.adapter_volume] * self.num_samples,
                        liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')
           
           # Add Ligation Master Mix (pre-mixed buffer + enzyme)
           pip_transfer(ham_int, self.tracked_tips_50uL, self.ligation_mastermix_positions, ligation_positions,
                        volumes=[self.ligation_mastermix_volume] * self.num_samples,
                        liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')
           
           # Mix samples
           self.tracked_tips_50uL.reset_all()  # Reset tips to ensure enough volume
           mix_plate(ham_int, self.tracked_tips_50uL, self.tip_support, self.num_samples,
                     self.HSP_Plate, mixing_volume=40, mix_cycles=10,
                     liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')
           
           # Incubate at room temperature for 15 minutes
           incubation_timer = start_timer(900)  # 15 minutes
           incubation_timer.wait(skip=self.device_simulation)

           print("Adapter Ligation completed.")

   def post_ligation_cleanup(self):
       """Step 4: Post-Ligation Cleanup with KAPA Pure Beads."""
       print("Starting Post-Ligation Cleanup step...")

       with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
           ham_int.initialize()
           normal_logging(ham_int, os.getcwd())

           # Transfer samples to MIDI plate
           current_volume = (self.dna_input_volume + self.fragmentation_mastermix_volume + 
                           self.end_repair_mastermix_volume + self.adapter_volume + 
                           self.ligation_mastermix_volume)
           
           transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                       self.HSP_Plate, self.MIDI_OffMagnet, volume=current_volume,
                       liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
           
           # Add KAPA Pure Beads (1.0X ratio)
           midi_positions = [(self.MIDI_OffMagnet, i) for i in range(self.num_samples)]
           multi_dispense(ham_int, self.tracked_tips_300uL, self.kapa_pure_beads_positions, midi_positions,
                          volumes=[self.post_ligation_bead_volume] * self.num_samples,
                          liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
           
           # Mix beads with samples
           mix_plate(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                     self.MIDI_OffMagnet, mixing_volume=80, mix_cycles=15,
                     liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
           
           # Incubate 5 minutes
           incubate_timer = start_timer(300)
           incubate_timer.wait(skip=self.device_simulation)
           
           # Move to magnet
           transport_resource(ham_int, self.MIDI_OffMagnet, self.MIDI_OnMagnet,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
           
           # Wait for beads to settle
           settle_timer = start_timer(120)
           settle_timer.wait(skip=self.device_simulation)
           
           # Remove supernatant
           transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                       self.MIDI_OnMagnet, self.MIDI_OffMagnet, 
                       volume=current_volume + self.post_ligation_bead_volume,
                       liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                       aspiration_height=0.3)
           
           # Ethanol washes (2x with 200µL)
           magnet_positions = [(self.MIDI_OnMagnet, i) for i in range(self.num_samples)]
           for wash in range(2):
               # Add ethanol
               transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                           self.ethanol_plate, self.MIDI_OnMagnet, volume=self.ethanol_wash_volume,
                           liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                           aspiration_height=0.3)
               
               # Wait 30 seconds
               wait_timer = start_timer(30)
               wait_timer.wait(skip=self.device_simulation)
               
               # Remove ethanol
               transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                           self.MIDI_OnMagnet, self.MIDI_OffMagnet, volume=self.ethanol_wash_volume,
                           liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                           aspiration_height=0.3)
           
           # Air dry for 2 minutes
           dry_timer = start_timer(120)
           dry_timer.wait(skip=self.device_simulation)
           
           # Move off magnet for elution
           transport_resource(ham_int, self.MIDI_OnMagnet.layout_name(), self.MIDI_OffMagnet.layout_name(),
                             resource_type=GrippedResource.MIDI, core_gripper=True)
           
           # Add nuclease-free water for elution
           offmagnet_positions = [(self.MIDI_OffMagnet, i) for i in range(self.num_samples)]
           pip_transfer(ham_int, self.tracked_tips_50uL, self.nuclease_free_water_positions, offmagnet_positions,
                          volumes=[self.post_ligation_elution_volume] * self.num_samples,
                          liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')
           
           # Mix for elution
           mix_plate(ham_int, self.tracked_tips_50uL, self.tip_support, self.num_samples,
                     self.MIDI_OffMagnet, mixing_volume=20, mix_cycles=15,
                     liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')
           
           # Incubate 2 minutes
           incubate_timer = start_timer(120)
           incubate_timer.wait(skip=self.device_simulation)

           # Move to magnet for final collection
           transport_resource(ham_int, self.MIDI_OffMagnet, self.MIDI_OnMagnet,
                             resource_type=GrippedResource.MIDI, core_gripper=True)
           
           # Wait for beads to settle
           settle_timer = start_timer(120)
           settle_timer.wait(skip=self.device_simulation)
           
           # Transfer eluted library to HSP_Plate_2 for PCR
           transfer_96(ham_int, self.tracked_tips_50uL, self.tip_support, self.num_samples,
                       self.MIDI_OnMagnet, self.HSP_Plate_2, volume=self.post_ligation_elution_volume - 2,
                       liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty',
                       aspiration_height=0.3)

           print("Post-Ligation Cleanup completed.")

   def library_amplification(self):
       """Step 5: Library Amplification with KAPA HiFi."""
       print("Starting Library Amplification step...")

       with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, persistent=self.persistent) as ham_int:
           ham_int.initialize()
           normal_logging(ham_int, os.getcwd())

           pcr_positions = [(self.HSP_Plate_2, i) for i in range(self.num_samples)]
           
           # Add KAPA HiFi HotStart ReadyMix (cold reagent)
           pip_transfer(ham_int, self.tracked_tips_50uL, self.kapa_hifi_mix_positions, pcr_positions,
                        volumes=[self.pcr_mix_volume] * self.num_samples,
                        liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')
           
           # Add Library Amplification Primer Mix (index primers)
           transfer_96(ham_int, self.tracked_tips_50uL, self.tip_support, self.num_samples,
                        self.HSP_Plate, self.HSP_Plate_2, volume=self.index_primer_volume,
                        liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')
           
           # Mix PCR reactions
           mix_plate(ham_int, self.tracked_tips_50uL, self.tip_support, self.num_samples,
                     self.HSP_Plate_2, mixing_volume=40, mix_cycles=5,
                     liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')
           
           # Move to thermal cycler for amplification
           odtc_open_door(ham_int, device_id=1)
           transport_resource(ham_int, self.HSP_Plate_2, self.HSP_ODTC,
                             resource_type=GrippedResource.PCR, core_gripper=True)
           
           # Add lid
           transport_resource(ham_int, self.Lid_Stack.fetch_next(), self.ODTC_Lid,
                             resource_type=GrippedResource.LID, core_gripper=True)
           
           odtc_close_door(ham_int, device_id=1)
           
           # Run PCR protocol (cycles based on DNA input)
           odtc_execute_protocol(ham_int, device_id=1, method_name='Library_Amplification_PCR.xml', simulating=self.device_simulation)

           odtc_wait_for_idle(ham_int, device_id=1, check_interval=30, max_wait=7200, simulating=self.device_simulation)
           
           # Remove from thermal cycler
           odtc_open_door(ham_int, device_id=1)
           transport_resource(ham_int, self.ODTC_Lid, self.Lid_Stack.put_back(),
                             resource_type=GrippedResource.LID, core_gripper=True)
           transport_resource(ham_int, self.HSP_ODTC, self.HSP_Plate_2,
                             resource_type=GrippedResource.PCR, core_gripper=True)
           odtc_close_door(ham_int, device_id=1)

           print("Library Amplification completed.")

   def final_cleanup_size_selection(self):
       """Step 6: Final Cleanup & Size Selection with KAPA Pure Beads."""
       print("Starting Final Cleanup & Size Selection step...")

       with HamiltonInterface(simulating=self.simulation, windowed=self.windowed, server_mode=False, persistent=self.persistent) as ham_int:
           ham_int.initialize()
           normal_logging(ham_int, os.getcwd())

           # Transfer PCR products to MIDI plate
           pcr_volume = (self.post_ligation_elution_volume - 2 + self.pcr_mix_volume + self.index_primer_volume)
           
           self.tracked_tips_300uL.reset_all()  # Reset tips to ensure enough volume
           transfer_96(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                       self.HSP_Plate_2, self.MIDI_OffMagnet, volume=pcr_volume,
                       liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')
           
           # Add KAPA Pure Beads (0.9X ratio for size selection)
           midi_positions = [(self.MIDI_OffMagnet, i) for i in range(self.num_samples)]
           pip_transfer(ham_int, self.tracked_tips_300uL, self.kapa_pure_beads_positions, midi_positions,
                          volumes=[self.final_cleanup_bead_volume] * self.num_samples,
                          liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')
           
           # Mix beads with PCR products
           mix_plate(ham_int, self.tracked_tips_300uL, self.tip_support, self.num_samples,
                     self.MIDI_OffMagnet, mixing_volume=80, mix_cycles=15, liquid_height=1,
                     liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
           
           # Incubate 5 minutes
           incubate_timer = start_timer(300)
           incubate_timer.wait(skip=self.device_simulation)
           
           # Move to magnet
           transport_resource(ham_int, self.MIDI_OffMagnet, self.MIDI_OnMagnet,
                             resource_type=GrippedResource.MIDI, core_gripper=True)


if __name__ == "__main__":
   protocol = KAPAHyperPlusProtocol(num_samples=96, dna_input_ng=100, fragmentation_time_min=15, pcr_cycles=8, 
                                    simulation=True, device_simulation=True)
   protocol.run_protocol()
   generate_tadm_report()