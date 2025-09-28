from pyhamilton import HamiltonInterface, LayoutManager, Reservoir60mL, Tip96, layout_item, Plate96
from pyhamilton_advanced import GRIPPER_CONFIGS, transport_resource, GripDirection, GrippedResource, GripperParams



lmgr = LayoutManager('PacBio_MultiPlexLibraryPrepDeck_v1.2.lay')

MIDI_CPAC01 = layout_item(lmgr, Plate96, 'MIDI_CPAC01')
MIDI_CPAC02 = layout_item(lmgr, Plate96, 'MIDI_CPAC02')

MIDI_Pipette_1 = layout_item(lmgr, Plate96, 'MIDI_Pipette01')
MIDI_Pipette_3 = layout_item(lmgr, Plate96, 'MIDI_Pipette03')

PCR_Pipette_1 = layout_item(lmgr, Plate96, 'PCR_Pipette01')
PCR_Pipette_2 = layout_item(lmgr, Plate96, 'PCR_Pipette02')
PCR_Pipette_3 = layout_item(lmgr, Plate96, 'PCR_Pipette03')


PCR_OnMagnet = layout_item(lmgr, Plate96, 'PCR_OnMagnet')
MIDI_OnMagnet = layout_item(lmgr, Plate96, 'MIDI_OnMagnet')

PCR_ODTC = layout_item(lmgr, Plate96, 'PCR_ODTC')
PCR_ODTC_inverted = layout_item(lmgr, Plate96, 'PCR_ODTC_inverted')

MIDI_Thermoshaker = layout_item(lmgr, Plate96, 'MIDI_ThermoShake')
MIDI_Thermoshaker_inverted = layout_item(lmgr, Plate96, 'MIDI_ThermoShake_inverted')

MIDI_Thermoshaker = layout_item(lmgr, Plate96, 'MIDI_ThermoShake_inverted')

with HamiltonInterface(windowed=True, simulating=False) as ham_int:
    ham_int.initialize()

    # MIDI Thermoshaker to MIDI on magnet
    transport_resource(ham_int, MIDI_Thermoshaker, MIDI_OnMagnet, iswap=True, grip_direction=GripDirection.LEFT, resource=GrippedResource.MIDI)
