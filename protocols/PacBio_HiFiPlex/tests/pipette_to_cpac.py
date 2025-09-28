#!/usr/bin/env python3
"""CPAC Pipetting Test"""

from pyhamilton import HamiltonInterface, LayoutManager, layout_item
from pyhamilton.resources import Plate96, TrackedTips, Tip96
from pyhamilton.pipetting import pip_transfer
from pyhamilton.consumables import ReagentTrackedEppiCarrier32

from pathlib import Path

def get_parent_lay_file():
    parent = Path(__file__).resolve().parent.parent  # go one level up
    for file in parent.glob("*.lay"):
        return str(file)  # return first match
    return None  # if no .lay file found


def test_cpac_pipetting(simulating=True):
    lmgr = LayoutManager(get_parent_lay_file())

    HSP_CPAC = layout_item(lmgr, Plate96, 'HSP_CPAC')
    
    # Assign small volume reagents to eppindorf tubes in carrier
    CAR_VIALS_SMALL = layout_item(lmgr, ReagentTrackedEppiCarrier32, 'CAR_VIALS_SMALL')
    FastSelect_position = CAR_VIALS_SMALL.assign_reagent_map('FastSelect', [11])
    RP_Primer_II_position = CAR_VIALS_SMALL.assign_reagent_map('RP_Primer_II', [1])
    
    # Tracked tips
    tracked_tips_50uL = TrackedTips.from_prefix(
        tracker_id="TIP_50ulF_L", volume_capacity=50, prefix="TIP_50ulF_L",
        count=8, tip_type=Tip96, lmgr=lmgr
    )
    
    num_samples = 8

    # Sample positions stored on CPAC plate
    HSP_CPAC_positions = [(HSP_CPAC, idx) for idx in range(num_samples)]
    volumes = [50] * num_samples
    
    windowed = not simulating
    with HamiltonInterface(simulating=simulating, windowed=windowed) as ham_int:
        ham_int.initialize()
        
        try:
            # Add FastSelect reagent from Eppindorf tube to CPAC sample positions
            pip_transfer(
                ham_int, tracked_tips_50uL, FastSelect_position,
                HSP_CPAC_positions, volumes,
                liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty',
                aspiration_height=1, dispense_height=1
            )

            # Add RP Primer II reagent from Eppindorf tube to CPAC sample positions
            pip_transfer(
                ham_int, tracked_tips_50uL, RP_Primer_II_position,
                HSP_CPAC_positions, volumes,
                liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty',
                aspiration_height=1, dispense_height=1
            )
            
            print("CPAC Pipetting: OK")
            
        except Exception as e:
            print(f"CPAC Pipetting: FAILED - {e}")

if __name__ == "__main__":
    test_cpac_pipetting()