from pyhamilton import (HamiltonInterface, LayoutManager, Reservoir60mL, TrackedTips, StackedResources, Tip96, Plate96, layout_item,
                        normal_logging)

import os

from pyhamilton_advanced import (shear_plate_96, double_aspirate_supernatant_96, ethanol_wash, pip_transfer, multi_dispense, 
                                    build_dispense_batches, batch_columnwise_positions, split_aspiration_positions)

lmgr = LayoutManager('PacBio_MultiPlexLibraryPrepDeck_v1.2.lay')

# Perhaps import stack management

MIDI_OffMagnet = layout_item(lmgr, Plate96, 'MIDI_Pipette')  # Assuming this is defined elsewhere in the layout
MagBeads_Container = layout_item(lmgr, Reservoir60mL, 'rgt_cont_60ml_BC_A00_0001')
ER_Mix = layout_item(lmgr, Reservoir60mL, 'rgt_cont_60ml_BC_A00_0002')
EDTA = layout_item(lmgr, Reservoir60mL, 'rgt_cont_60ml_BC_A00_0003')

tips = tip_tracker_50uL = TrackedTips.from_prefix(
    tracker_id="HTF_L",
    prefix="HTF_L",
    count=3,
    tip_type=Tip96,
    lmgr=lmgr,
    reset=True  # Reset the tracker state
)


def condense_volumes(lst, max_volume):
    total = sum(lst)
    return [max_volume] * (total // max_volume) + ([total % max_volume] if 0 < total % max_volume >= min(lst) else [])


dispense_positions = [(MIDI_OffMagnet, idx) for idx in range(96)]
dispense_volumes = [50]*96  # Assuming 50 uL for each sample

aspirate_positions = [(MagBeads_Container, idx) for idx in range(8)]

aspirate_volumes = condense_volumes(dispense_volumes, 200)
column_dispense_positions = batch_columnwise_positions(dispense_positions) # Batch dispense positions into columns of 8
column_volumes_list = batch_columnwise_positions(dispense_volumes)

column_aspiration_positions = batch_columnwise_positions(aspirate_positions) # Batch aspiration positions into columns of 8
column_aspiration_volumes_list = batch_columnwise_positions(aspirate_volumes)

column_aspiration_volumes = [300]*8
batched_vols = build_dispense_batches(column_aspiration_volumes, column_dispense_positions, column_volumes_list)
#print("Batched Volumes:", batched_vols)

with HamiltonInterface(windowed=True, simulating=False) as ham_int:
    ham_int.initialize()
    normal_logging(ham_int, os.getcwd())
    aspiration_positions = [(MagBeads_Container, idx) for idx in range(8)]
    dispense_positions = [(MIDI_OffMagnet, idx) for idx in range(96)]
    volumes = [50]*96
#
    multi_dispense(ham_int, tips, aspiration_positions, dispense_positions, volumes, 
                 liquid_class = 'HighVolumeFilter_Water_DispenseSurface_Empty')
    
    