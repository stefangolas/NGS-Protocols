from pyhamilton import (HamiltonInterface, LayoutManager, Tip96, TrackedTips, resource_list_with_prefix, 
                        Tip96, Plate96, ResourceType, tip_pick_up, DeckResource, 
                        tracked_tip_pick_up, tracked_tip_pick_up_96, StackedResources, move_plate_using_gripper)


BioRadHardShell_Stack = StackedResources.from_prefix("BioRadHardShell_Stack", "BioRadHardShell_Stack", 3)
AbgeneMIDI_Stack = StackedResources.from_prefix("AbgeneMIDI_Stack1", "AbgeneMIDI_Stack1", 4)


with HamiltonInterface(windowed=True) as ham_int:
    ham_int.initialize()
    for _ in range(3):
        plate_seq = AbgeneMIDI_Stack.fetch_next()
        move_plate_using_gripper(ham_int, plate_seq, 'MIDI_OnMagnet', gripHeight=5)
