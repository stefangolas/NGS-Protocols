"""10X GEX Library Prep - Standalone Test Script"""

import os
from pathlib import Path
from pyhamilton import HamiltonInterface, LayoutManager, layout_item, start_timer, normal_logging
from pyhamilton.pipetting import pip_transfer, transfer_96, pip_mix, mix_plate
from pyhamilton.consumables import (ReagentTrackedReservoir60mL, ReagentTrackedBulkPlate,
                                   ReagentTrackedEppiCarrier32, ReagentTrackedPlate96)
from pyhamilton.resources import Plate96, Tip96, Lid, Waste96, TrackedTips, TipSupportTracker, StackedResources
from pyhamilton.transport import transport_resource, GrippedResource
from pyhamilton.devices import (hhs_set_simulation, hhs_create_usb_device,
                            hhs_start_shaker, hhs_stop_shaker)


def get_parent_lay_file():
    """Find the .lay file in the parent directory."""
    parent = Path(__file__).resolve().parent.parent
    for file in parent.glob("*.lay"):
        return str(file)
    return None


def deck_setup(simulating=True):
    """
    Args:
        simulating (bool): Run Hamilton in simulation mode
        device_simulation (bool): Simulate HHS devices
    """
    # Load layout
    lay_file = get_parent_lay_file()
    lmgr = LayoutManager(lay_file)
    
    # Protocol parameters
    num_samples = 8
    
    # PCR plates for thermal cycling steps
    HSP_Plate = layout_item(lmgr, Plate96, 'HSP_Pipette')
    HSP_Plate2 = layout_item(lmgr, Plate96, 'HSP_Pipette2')
    HSP_Waste = layout_item(lmgr, Plate96, 'HSP_Waste')
    HHS1_HSP = layout_item(lmgr, Plate96, 'HHS1_HSP')
    HSP_Magnet = layout_item(lmgr, Plate96, 'HSP_OnMagnet')
    HSP_ODTC = layout_item(lmgr, Plate96, 'HSP_ODTC')
    HSP_ODTC_Lid = layout_item(lmgr, Lid, 'Ham_ComfortLid_ODTC')
    
    # MIDI plates for magnetic bead operations
    MIDI_Plate = layout_item(lmgr, Plate96, 'MIDI_Pipette')
    MIDI_Waste = layout_item(lmgr, Plate96, 'MIDI_Waste')
    MIDI_CPAC = layout_item(lmgr, Plate96, 'MIDI_CPAC')
    MIDI_OnMagnet = layout_item(lmgr, Plate96, 'MIDI_OnMagnet')
    HHS5_MIDI = layout_item(lmgr, Plate96, 'HHS5_MIDI')
    
    # Reagent containers
    CAR_VIALS_SMALL = layout_item(lmgr, ReagentTrackedEppiCarrier32, 'CAR_VIALS_SMALL')
    CPAC_Reagents = layout_item(lmgr, ReagentTrackedPlate96, 'CPAC_HSP_0001')
    RGT_01 = layout_item(lmgr, ReagentTrackedReservoir60mL, 'rgt_cont_60ml_BC_A00_0001')
    RGT_02 = layout_item(lmgr, ReagentTrackedReservoir60mL, 'rgt_cont_60ml_BC_A00_0002')
    EthanolReservoir = layout_item(lmgr, ReagentTrackedBulkPlate, 'RGT_Ethanol')
    
    # Reagent mapping and positions
    # Room temperature small volume reagents
    fragmentation_buffer_positions = CAR_VIALS_SMALL.assign_reagent_map('FragmentationBuffer', [0])
    ligation_mix_positions = CAR_VIALS_SMALL.assign_reagent_map('LigationMix', [1])
    buffer_eb_positions = CAR_VIALS_SMALL.assign_reagent_map('BufferEB', [2])
    
    # Cold-sensitive reagents in CPAC
    fragmentation_enzyme_positions = CPAC_Reagents.assign_reagent_map('FragmentationEnzyme', [0])
    dna_ligase_positions = CPAC_Reagents.assign_reagent_map('DNALigase', [1])
    library_amp_mix_positions = CPAC_Reagents.assign_reagent_map('LibraryAmpMix', [2])
    

    # Bulk reagents
    nuclease_free_water_positions = RGT_01.assign_reagent_map('NucleaseFreeWater', range(8))
    spriselect_positions = RGT_02.assign_reagent_map('SPRIselect', range(8))
    
    # Index positions
    index_positions = [(HHS1_HSP, i) for i in range(num_samples)]
    
    # Liquid waste
    Liquid_Waste = layout_item(lmgr, Waste96, 'core96externalwaste_0001')
    
    # Stacked resources
    HSP_Stack = StackedResources.from_prefix(
        tracker_id="BioRadHardShell_Stack4",
        prefix="BioRadHardShell_Stack4",
        lmgr=lmgr,
        resource_type=Plate96,
        count=5)
    
    Lid_Stack = StackedResources.from_prefix(
        tracker_id="Ham_ComfortLid_Stack",
        prefix="Ham_ComfortLid_Stack",
        lmgr=lmgr,
        resource_type=Lid,
        count=3)
    
    MIDI_Stack = StackedResources.from_prefix(
        tracker_id="AbgeneMIDI_Stack1",
        prefix="AbgeneMIDI_Stack1",
        lmgr=lmgr,
        resource_type=Plate96,
        count=3)
    
    # Tracked tips
    tracked_tips_50uL = TrackedTips.from_prefix(
        tracker_id="TIP_50uLF_L",
        volume_capacity=50,
        prefix="TIP_50uLF_L",
        count=8,
        tip_type=Tip96,
        lmgr=lmgr)
    
    tracked_tips_300uL = TrackedTips.from_prefix(
        tracker_id="STF_L",
        volume_capacity=300,
        prefix="STF_L",
        count=8,
        tip_type=Tip96,
        lmgr=lmgr)
    
    tracked_tips_1000uL = TrackedTips.from_prefix(
        tracker_id="HTF_L",
        volume_capacity=1000,
        prefix="HTF_L",
        count=2,
        tip_type=Tip96,
        lmgr=lmgr)
    
    # Tip support
    tip_support_resource = layout_item(lmgr, Tip96, 'TipSupport_0001')
    tip_support = TipSupportTracker(tip_support_resource)
    
    # Main protocol execution
    windowed = not simulating
    with HamiltonInterface(simulating=simulating, windowed=windowed) as ham_int:
        ham_int.initialize()
        normal_logging(ham_int, os.getcwd())
        
        print("Deck layout initialized successfully")


if __name__ == "__main__":
    # Run with simulation mode
    # Set simulating=False to run on real hardware
    # Set device_simulation=False to use real HHS devices
    deck_setup(simulating=True, device_simulation=True)