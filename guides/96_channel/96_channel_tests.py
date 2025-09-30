import os
from pathlib import Path
from pyhamilton import (HamiltonInterface, LayoutManager, layout_item, 
                        TrackedTips, TipSupportTracker, normal_logging)
from pyhamilton.pipetting import transfer_96
from pyhamilton.resources import Plate96, Tip96


def get_parent_lay_file():
    """Find the .lay file in the parent directory."""
    parent = Path(__file__).resolve().parent.parent
    for file in parent.glob("*.lay"):
        return str(file)
    return None


def mph_tip_support(simulating=True):
    """
    Demonstrate automatic tip volume switching.
    Three transfers: 300 µL tips → 50 µL tips → 300 µL tips
    """
    
    # Load layout
    lay_file = get_parent_lay_file()
    lmgr = LayoutManager(lay_file)
    
    # Setup plates
    source_plate = layout_item(lmgr, Plate96, 'HSP_Pipette')
    dest_plate = layout_item(lmgr, Plate96, 'HSP_Pipette2')
    
    # Setup 300 µL tips
    tracked_tips_300uL = TrackedTips.from_prefix(
        tracker_id="stf_l",
        volume_capacity=300,
        prefix="stf_l",
        count=8,
        tip_type=Tip96,
        lmgr=lmgr
    )
    
    # Setup 50 µL tips
    tracked_tips_50uL = TrackedTips.from_prefix(
        tracker_id="TIP_50ulF_L",
        volume_capacity=50,
        prefix="TIP_50ulF_L",
        count=8,
        tip_type=Tip96,
        lmgr=lmgr
    )
    
    # Setup tip support
    tip_support_resource = layout_item(lmgr, Tip96, 'TipSupport_0001')
    tip_support = TipSupportTracker(tip_support_resource)
    
    # Initialize Hamilton
    windowed = not simulating
    with HamiltonInterface(simulating=simulating, windowed=windowed) as ham_int:
        ham_int.initialize()
        normal_logging(ham_int, os.getcwd())
        
        
        transfer_96(
            ham_int, tracked_tips_300uL, tip_support,
            num_samples=8,
            source_plate=source_plate,
            target_plate=dest_plate,
            volume=100,
            liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty'
        )
        
        transfer_96(
            ham_int, tracked_tips_50uL, tip_support,
            num_samples=8,
            source_plate=source_plate,
            target_plate=dest_plate,
            volume=20,
            liquid_class='Tip_50ulFilter_Water_DispenseJet_Empty'
        )
        transfer_96(
            ham_int, tracked_tips_300uL, tip_support,
            num_samples=8,
            source_plate=source_plate,
            target_plate=dest_plate,
            volume=150,
            liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty'
        )
        print("  ✓ Complete\n")
        
        print("=== Test Complete ===")
        print("TipSupport automatically swapped racks when tip volume changed")


if __name__ == "__main__":
    mph_tip_support(simulating=True)