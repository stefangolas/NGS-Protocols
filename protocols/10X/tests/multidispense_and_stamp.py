from pyhamilton import HamiltonInterface, LayoutManager, layout_item
from pyhamilton.pipetting import pip_mix, multi_dispense, transfer_96
from pyhamilton.consumables import ReagentTrackedEppiCarrier32, ReagentTrackedPlate96
from pyhamilton.resources import Plate96, Tip96, TrackedTips, TipSupportTracker
from pyhamilton.transport import transport_resource, GrippedResource
from pathlib import Path


def get_parent_lay_file():
    parent = Path(__file__).resolve().parent.parent
    for file in parent.glob("*.lay"):
        return str(file)
    return None


def test_multidispense_and_transfer(simulating=True):
    lmgr = LayoutManager(get_parent_lay_file())

    # Setup resources exactly as in protocol
    CAR_VIALS_SMALL = layout_item(lmgr, ReagentTrackedEppiCarrier32, 'CAR_VIALS_SMALL')
    CPAC_Reagents = layout_item(lmgr, ReagentTrackedPlate96, 'CPAC_HSP_0001')
    MIDI_CPAC = layout_item(lmgr, Plate96, 'MIDI_CPAC')
    MIDI_Plate = layout_item(lmgr, Plate96, 'MIDI_Plate')
    HSP_Plate2 = layout_item(lmgr, Plate96, 'HSP_Plate2')

    # Reagent mapping
    fragmentation_buffer_positions = CAR_VIALS_SMALL.assign_reagent_map('FragmentationBuffer', [0])
    buffer_eb_positions = CAR_VIALS_SMALL.assign_reagent_map('BufferEB', [2])
    fragmentation_enzyme_positions = CPAC_Reagents.assign_reagent_map('FragmentationEnzyme', [0])

    # Tips
    tracked_tips_50uL = TrackedTips.from_prefix(
        tracker_id="TIP_50uLF_L", volume_capacity=50, prefix="TIP_50uLF_L",
        count=8, tip_type=Tip96, lmgr=lmgr
    )
    tracked_tips_300uL = TrackedTips.from_prefix(
        tracker_id="STF_L", volume_capacity=300, prefix="STF_L",
        count=8, tip_type=Tip96, lmgr=lmgr
    )
    tip_support = TipSupportTracker(layout_item(lmgr, Tip96, 'TipSupport_0001'))

    # Protocol parameters
    num_samples = 8
    excess_factor = 1.1
    total_reactions = int(num_samples * excess_factor)

    # Volume calculations
    buffer_eb_vol = 20 * total_reactions
    frag_buffer_vol = 5 * total_reactions
    frag_enzyme_vol = 10 * total_reactions

    windowed = not simulating
    with HamiltonInterface(simulating=simulating, windowed=windowed) as ham_int:
        ham_int.initialize()

        try:
            # Step 1: Mix fragmentation buffer first (room temp reagent)
            pip_mix(ham_int, tracked_tips_50uL, fragmentation_buffer_positions,
                   mix_volume=200, mix_cycles=10,
                   liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')

            # Step 2: Setup mix positions in MIDI plate on CPAC
            mix_position = [(MIDI_CPAC, idx) for idx in range(num_samples)]

            # Step 3: Add Buffer EB (room temp)
            multi_dispense(ham_int, tips=tracked_tips_300uL, source_positions=buffer_eb_positions,
                          dispense_positions=mix_position, volumes=[buffer_eb_vol] * num_samples,
                          aspiration_height=0, liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')

            # Step 4: Add Fragmentation Buffer (room temp)
            multi_dispense(ham_int, tips=tracked_tips_300uL, source_positions=fragmentation_buffer_positions,
                          dispense_positions=mix_position, volumes=[frag_buffer_vol] * num_samples,
                          aspiration_height=0, liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')

            # Step 5: Add Fragmentation Enzyme (cold reagent from CPAC)
            multi_dispense(ham_int, tracked_tips_300uL, fragmentation_enzyme_positions, mix_position,
                          volumes=[frag_enzyme_vol] * num_samples, aspiration_height=0,
                          liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')

            # Step 6: Mix the prepared fragmentation mix
            total_vol = buffer_eb_vol + frag_buffer_vol + frag_enzyme_vol
            pip_mix(ham_int, tips=tracked_tips_300uL, positions_to_mix=mix_position,
                   mix_cycles=3, mix_volume=25, liquid_height=0,
                   liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')

            # Step 7: Transport MIDI plate from CPAC to MIDI_Plate
            transport_resource(ham_int, MIDI_CPAC, MIDI_Plate,
                             resource_type=GrippedResource.MIDI, core_gripper=True)

            # Step 8: Transfer fragmentation mix to HSP_Plate2
            transfer_96(ham_int, tracked_tips_300uL, tip_support, num_samples,
                       source_plate=MIDI_Plate, target_plate=HSP_Plate2, volume=40,
                       liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')

            print("Multidispense and Transfer: OK")

        except Exception as e:
            print(f"Multidispense and Transfer: FAILED - {e}")


if __name__ == "__main__":
    test_multidispense_and_transfer(simulating=True)