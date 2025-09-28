#!/usr/bin/env python3
"""Magnetic Bead Cleanup Test - Exact Protocol Sequence"""

from pyhamilton import HamiltonInterface, LayoutManager, layout_item, start_timer
from pyhamilton.pipetting import pip_transfer, transfer_96, mix_plate, double_aspirate_supernatant_96
from pyhamilton.consumables import ReagentTrackedReservoir60mL, ReagentTrackedBulkPlate
from pyhamilton.resources import Plate96, Tip96, TrackedTips, TipSupportTracker, StackedResources
from pyhamilton.transport import transport_resource, GrippedResource
from pathlib import Path


def get_parent_lay_file():
    parent = Path(__file__).resolve().parent.parent
    for file in parent.glob("*.lay"):
        return str(file)
    return None


def test_magnetic_bead_cleanup(simulating=True):
    lmgr = LayoutManager(get_parent_lay_file())

    # Setup resources exactly as in protocol
    HSP_Pipette = layout_item(lmgr, Plate96, 'HSP_Pipette')
    HSP_Pipette2 = layout_item(lmgr, Plate96, 'HSP_Pipette2')
    HSP_Waste = layout_item(lmgr, Plate96, 'HSP_Waste')
    MIDI_Pipette = layout_item(lmgr, Plate96, 'MIDI_Pipette')
    MIDI_OnMagnet = layout_item(lmgr, Plate96, 'MIDI_OnMagnet')
    MPH_Waste = layout_item(lmgr, Plate96, 'MPH_Waste')
    
    HSP_Stack = StackedResources.from_prefix(
        tracker_id="BioRadHardshell_Stack1", prefix="BioRadHardshell_Stack1",
        count=4, lmgr=lmgr, resource_type=Plate96
    )
    
    RGT_01 = layout_item(lmgr, ReagentTrackedReservoir60mL, 'RGT_01')
    EthanolReservoir = layout_item(lmgr, ReagentTrackedBulkPlate, 'RGT_Ethanol')
    ampure_bead_positions = RGT_01.assign_reagent_map('AmpureBead', range(8))
    elution_buffer_positions = RGT_01.assign_reagent_map('ElutionBuffer', range(8))
    
    tracked_tips_300uL = TrackedTips.from_prefix(
        tracker_id="stf_l", volume_capacity=300, prefix="stf_l",
        count=8, tip_type=Tip96, lmgr=lmgr
    )
    tip_support = TipSupportTracker(layout_item(lmgr, Tip96, 'TipSupport_0001'))
    
    # Protocol parameters from original
    num_samples = 8
    magbead_volume = 50
    elution_volume = 30
    adapter_ligation_mix_volume = 25
    adapter_mix_volume = 15
    total_ligation_volume = elution_volume - 2 + adapter_ligation_mix_volume + adapter_mix_volume
    supernatant_removal_volume = 150
    fragment_buffer_volume = 200
    
    windowed = not simulating
    with HamiltonInterface(simulating=simulating, windowed=windowed) as ham_int:
        ham_int.initialize()
        
        try:
            # Step 1: Add magnetic beads to MIDI positions
            midi_positions = [(MIDI_Pipette, i) for i in range(num_samples)]
            pip_transfer(ham_int, tracked_tips_300uL, ampure_bead_positions, midi_positions,
                        volumes=[magbead_volume] * num_samples,
                        liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                        aspiration_height=1, dispense_height=1)

            # Step 2: Transfer ligation reaction to MIDI plate
            transfer_96(ham_int, tips=tracked_tips_300uL, tip_support=tip_support, num_samples=num_samples,
                       source_plate=HSP_Pipette, target_plate=MIDI_Pipette,
                       volume=total_ligation_volume,
                       liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                       aspiration_height=1, dispense_height=1)

            # Step 3: Move HSP plate to waste
            transport_resource(ham_int, HSP_Pipette, HSP_Waste, 
                             core_gripper=True, resource_type=GrippedResource.PCR)

            # Step 4: Mix and incubate
            mix_plate(ham_int, tracked_tips_300uL, tip_support, num_samples,
                     MIDI_Pipette, mixing_volume=50, mix_cycles=5,
                     liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            incubation_timer = start_timer(300)
            incubation_timer.wait(skip=simulating)

            # Step 5: Move to magnet
            transport_resource(ham_int, MIDI_Pipette, MIDI_OnMagnet, 
                             core_gripper=True, resource_type=GrippedResource.PCR)

            settle_timer = start_timer(120)
            settle_timer.wait(skip=simulating)

            # Step 6: Remove supernatant
            double_aspirate_supernatant_96(ham_int, tracked_tips_300uL, tip_support, num_samples,
                                          MIDI_OnMagnet, MPH_Waste,
                                          first_volume=supernatant_removal_volume,
                                          second_volume=30,
                                          liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                                          first_aspiration_height=0.5)

            # Step 7: 2x Fragment buffer washes
            for wash_cycle in range(2):
                # Add fragment buffer
                transfer_96(ham_int, tracked_tips_300uL, tip_support, num_samples,
                           EthanolReservoir, MIDI_OnMagnet,
                           volume=fragment_buffer_volume,
                           liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                           aspiration_height=1, dispense_height=1)

                # Mix
                mix_plate(ham_int, tracked_tips_300uL, tip_support, num_samples,
                          MIDI_OnMagnet, mixing_volume=40, mix_cycles=5,
                          liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

                # Incubate
                wash_timer = start_timer(120)
                wash_timer.wait(skip=simulating)

                # Remove supernatant
                double_aspirate_supernatant_96(ham_int, tracked_tips_300uL, tip_support, num_samples,
                                              MIDI_OnMagnet, MPH_Waste,
                                              first_volume=fragment_buffer_volume + 20,
                                              second_volume=30,
                                              liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                                              first_aspiration_height=0.5)

            # Step 8: Air dry beads
            dry_timer = start_timer(120)
            dry_timer.wait(skip=simulating)

            # Step 9: Get fresh HSP plate for final elution
            transport_resource(ham_int, HSP_Stack.fetch_next(), HSP_Pipette2, 
                             core_gripper=True, resource_type=GrippedResource.PCR)

            # Step 10: Final elution with Elution Buffer
            MIDI_OnMagnet_positions = [(MIDI_OnMagnet, i) for i in range(num_samples)]
            pip_transfer(ham_int, tracked_tips_300uL,
                       elution_buffer_positions, MIDI_OnMagnet_positions,
                       volumes=[elution_volume] * num_samples,
                       liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                       aspiration_height=1, dispense_height=1)

            # Step 11: Mix and incubate
            mix_plate(ham_int, tracked_tips_300uL, tip_support, num_samples,
                     MIDI_OnMagnet, mixing_volume=12, mix_cycles=5,
                     liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            final_mix_timer = start_timer(120)
            final_mix_timer.wait(skip=simulating)

            # Step 12: Collect final eluted DNA
            transfer_96(ham_int, tracked_tips_300uL, tip_support, num_samples,
                       MIDI_OnMagnet, HSP_Pipette2,
                       volume=elution_volume - 2,  # Leave 2µL to avoid beads
                       liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                       aspiration_height=0.5, dispense_height=1)
            
            print("Magnetic Bead Cleanup: OK")
            
        except Exception as e:
            print(f"Magnetic Bead Cleanup: FAILED - {e}")


if __name__ == "__main__":
    test_magnetic_bead_cleanup(simulating=True)