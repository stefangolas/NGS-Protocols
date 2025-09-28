from pyhamilton import (HamiltonInterface, LayoutManager, Tip96, TrackedTips, resource_list_with_prefix,
                        Tip96, Plate96, ResourceType,    StackedResources, layout_item)

from pyhamilton.transport import transport_resource, GrippedResource


def get_parent_lay_file():
    from pathlib import Path
    parent = Path(__file__).resolve().parent.parent  # go one level up
    for file in parent.glob("*.lay"):
        return str(file)  # return first match
    return None  # if no .lay file found

lmgr = LayoutManager(get_parent_lay_file())

HSP_Stack = StackedResources.from_prefix(
    tracker_id="BioRadHardshell_Stack1",
    prefix="BioRadHardshell_Stack1",
    count=4,
    lmgr=lmgr,
    resource_type=Plate96
    )

MIDI_Stack = StackedResources.from_prefix(
    tracker_id="ABGENE_MIDI_Stack1",
    prefix="ABGENE_MIDI_Stack1",
    count=4,
    lmgr=lmgr,
    resource_type=Plate96
)


HSP_Pipette = layout_item(lmgr, Plate96, 'HSP_Pipette')
HHS3_MIDI = layout_item(lmgr, Plate96, 'HHS3_MIDI')

with HamiltonInterface(windowed=True) as ham_int:
    ham_int.initialize()
    transport_resource(ham_int, HSP_Stack.fetch_next(), HSP_Pipette, resource_type=GrippedResource.PCR, core_gripper=True)
    transport_resource(ham_int, MIDI_Stack.fetch_next(), HHS3_MIDI, resource_type=GrippedResource.PCR, core_gripper=True)