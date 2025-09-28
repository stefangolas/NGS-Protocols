from pyhamilton import (HamiltonInterface, LayoutManager, Reservoir60mL, TrackedTips, StackedResources, Tip96, Plate96, layout_item,
                        normal_logging, hhs_create_star_device, hhs_set_temp_param, hhs_start_temp_ctrl, hhs_stop_temp_ctrl, hhs_set_simulation)

import os

from pyhamilton_advanced import shear_plate_96, double_aspirate_supernatant_96, ethanol_wash, pip_transfer, mix_plate

lmgr = LayoutManager('PacBio_MultiPlexLibraryPrepDeck_v1.2.lay')

# Perhaps import stack management

MIDI_OffMagnet = layout_item(lmgr, Plate96, 'MIDI_Pipette')  # Assuming this is defined elsewhere in the layout
MagBeads_Container = layout_item(lmgr, Reservoir60mL, 'rgt_cont_60ml_BC_A00_0001')
ER_Mix = layout_item(lmgr, Reservoir60mL, 'rgt_cont_60ml_BC_A00_0002')
EDTA = layout_item(lmgr, Reservoir60mL, 'rgt_cont_60ml_BC_A00_0003')

tips = tip_tracker_50uL = TrackedTips.from_prefix(
    tracker_id="TIP_50uLF_L",
    prefix="TIP_50uLF_L",
    count=8,
    tip_type=Tip96,
    lmgr=lmgr,
    reset=True  # Reset the tracker state
)


# This works
with HamiltonInterface(windowed=True, simulating=False) as ham_int:
    ham_int.initialize()
    normal_logging(ham_int, os.getcwd())
    
    hhs_set_simulation(ham_int, 1)  # Set simulation mode if needed
    hhs_create_star_device(ham_int, 'ML_STAR', 1)
    hhs_start_temp_ctrl(ham_int, 1, 37, True)
    hhs_stop_temp_ctrl(ham_int, 1)