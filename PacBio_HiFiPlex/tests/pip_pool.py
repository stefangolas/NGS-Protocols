from pyhamilton import (HamiltonInterface, LayoutManager, Reservoir60mL, TrackedTips, StackedResources, Tip96, Plate96, layout_item,
                        normal_logging)
from pyhamilton.pipetting import pip_pool
from pyhamilton.consumables import ReagentTrackedFalconCarrier24
import os


lmgr = LayoutManager('PacBio_MultiPlexLibraryPrepDeck_v1.2.lay')

# Perhaps import stack management

MIDI_OffMagnet = layout_item(lmgr, Plate96, 'MIDI_Pipette')  # Assuming this is defined elsewhere in the layout
PoolingTubes = layout_item(lmgr, ReagentTrackedFalconCarrier24, 'SMP_CAR_24_15x75_A00_0001')
pooling_positions = [(PoolingTubes, idx) for idx in range(1)]

tips = tip_tracker_50uL = TrackedTips.from_prefix(
    tracker_id="TIP_50uLF_L",
    prefix="TIP_50uLF_L",
    volume_capacity=50,
    count=8,
    tip_type=Tip96,
    lmgr=lmgr,
    reset=True  # Reset the tracker state
)


# This seems to work. We don't yet query the user for refilling troughs or try to accumulate residual volumes.
with HamiltonInterface(windowed=False, simulating=True) as ham_int:
    ham_int.initialize()
    normal_logging(ham_int, os.getcwd())

    aspiration_positions = [(MIDI_OffMagnet, idx) for idx in range(96)]
    volumes = [50]*96

    pip_pool(ham_int, tips, aspiration_positions, pooling_positions, volumes, 
                 liquid_class = 'Tip_50ulFilter_Water_DispenseSurface_Empty', aspiration_height=1,
                 dispense_height=1)