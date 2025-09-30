"""Sample Cleanup 1 Test - Complete Protocol with HHS Support"""

import os
from pyhamilton import HamiltonInterface, LayoutManager, layout_item, start_timer, normal_logging
from pyhamilton.pipetting import pip_transfer, transfer_96, pip_mix, double_aspirate_supernatant_96, ethanol_wash
from pyhamilton.consumables import ReagentTrackedReservoir60mL, ReagentTrackedBulkPlate
from pyhamilton.resources import Plate96, Tip96, TrackedTips, TipSupportTracker, StackedResources
from pyhamilton.transport import transport_resource, GrippedResource
from pyhamilton.devices import (hhs_set_simulation, hhs_create_usb_device, 
                            hhs_start_shaker, hhs_stop_shaker)
from pathlib import Path


def get_parent_lay_file():
    """Find the .lay file in the parent directory."""
    parent = Path(__file__).resolve().parent.parent
    for file in parent.glob("*.lay"):
        return str(file)
    return None


class HHS:
    """Hamilton Heater Shaker wrapper class."""
    def __init__(self, node, sequence, lmgr):
        self.node = node
        self.sequence = sequence
        self.resource = layout_item(lmgr, Plate96, sequence)



def magnetic_bead_cleanup(simulating=True, device_simulation=True):
    """
    Execute Sample Cleanup 1 protocol with two rounds of magnetic bead purification.
    
    Args:
        simulating (bool): Run Hamilton in simulation mode
        device_simulation (bool): Simulate HHS devices (timers still run)
    """
    # Load layout
    lay_file = get_parent_lay_file()
    lmgr = LayoutManager(lay_file)
    
    # Protocol parameters
    num_samples = 8
    sample_volume = 50
    post_shear_magbead_volume = sample_volume
    post_shear_etoh_wash_volume = 200
    post_shear_elution_buffer_volume = 30
    post_shear_elution_volume = 25.5
    
    # Setup HSP plates
    HSP_CPAC = layout_item(lmgr, Plate96, 'HSP_CPAC')
    HSP_Pipette2 = layout_item(lmgr, Plate96, 'HSP_Pipette2')
    HSP_Waste = layout_item(lmgr, Plate96, 'HSP_Waste')
    
    # Setup MIDI plates
    MIDI_OnMagnet = layout_item(lmgr, Plate96, 'MIDI_OnMagnet')
    MIDI_Waste = layout_item(lmgr, Plate96, 'MIDI_Waste')

    # MPH Waste
    MPH_Waste = layout_item(lmgr, Plate96, 'MPH_Waste')
    
    # Setup HHS devices
    HHS3_MIDI = HHS(node=3, sequence="HHS3_MIDI", lmgr=lmgr)
    HHS5_MIDI = HHS(node=5, sequence="HHS5_MIDI", lmgr=lmgr)
    
    # Setup reagents
    RGT_01 = layout_item(lmgr, ReagentTrackedReservoir60mL, 'RGT_01')
    EthanolReservoir = layout_item(lmgr, ReagentTrackedBulkPlate, 'RGT_Ethanol')
    
    QIAseq_Beads_positions = RGT_01.assign_reagent_map('QIAseq_Beads', range(8))
    Nuclease_Free_Water_positions = RGT_01.assign_reagent_map('Nuclease_Free_Water', range(8))
    
    # Setup stacked resources
    HSP_Stack = StackedResources.from_prefix(
        tracker_id="BioRadHardshell_Stack1", prefix="BioRadHardshell_Stack1",
        count=4, lmgr=lmgr, resource_type=Plate96
    )
    
    MIDI_Stack = StackedResources.from_prefix(
        tracker_id="ABGENE_MIDI_Stack1", prefix="ABGENE_MIDI_Stack1",
        count=4, lmgr=lmgr, resource_type=Plate96
    )
    
    # Setup tips
    tracked_tips_300uL = TrackedTips.from_prefix(
        tracker_id="stf_l", volume_capacity=300, prefix="stf_l",
        count=8, tip_type=Tip96, lmgr=lmgr
    )
    
    tip_support = TipSupportTracker(layout_item(lmgr, Tip96, 'TipSupport_0001'))
    
    # Main protocol execution
    windowed = not simulating
    with HamiltonInterface(simulating=simulating, windowed=windowed) as ham_int:
        ham_int.initialize()
        normal_logging(ham_int, os.getcwd())
        
        # Initialize HHS devices
        print("Initializing HHS devices...")
        hhs_set_simulation(ham_int, 1 if device_simulation else 0)
        for node in [3, 5]:
            try:
                hhs_create_usb_device(ham_int, node)
                print(f"  HHS node {node}: OK")
            except Exception as e:
                print(f"  Warning: Could not initialize HHS node {node}: {e}")
        
        print("\n=== Starting Sample Cleanup 1 Protocol ===\n")
        
        # Get fresh MIDI plate from stack and move to HHS3
        print("Step 1: Preparing MIDI plate on HHS3...")
        transport_resource(ham_int, MIDI_Stack.fetch_next(), HHS3_MIDI.resource,
                         resource_type=GrippedResource.MIDI, core_gripper=True)
        
        # Add QIAseq Beads to HHS3 position (Round 1)
        print("Step 2: Adding QIAseq Beads (Round 1)...")
        HHS3_MIDI_positions = [(HHS3_MIDI.resource, idx) for idx in range(num_samples)]
        volumes = [post_shear_magbead_volume] * num_samples

        # Mix beads thoroughly before transfering to MIDI plate
        pip_mix(ham_int, tips=tracked_tips_300uL, positions_to_mix=QIAseq_Beads_positions,
                mix_volume=post_shear_elution_volume, mix_cycles=10,
                liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                liquid_height=1)

        pip_transfer(ham_int, tracked_tips_300uL, QIAseq_Beads_positions,
                    HHS3_MIDI_positions, volumes,
                    liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                    aspiration_height=1, dispense_height=1)
        
        # Transfer samples from HSP_CPAC to MIDI plate
        print("Step 3: Transferring samples to MIDI plate...")
        transfer_96(ham_int, tracked_tips_300uL, tip_support=tip_support, num_samples=num_samples,
                   target_plate=HSP_CPAC, source_plate=HHS3_MIDI.resource, volume=sample_volume,
                   aspiration_mix_cycles=3, aspiration_mix_volume=sample_volume,
                   liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                   aspiration_height=1, dispense_height=1)
        
        # Shake plate
        print("Step 4: Mixing on shaker (10 sec @ 1000 rpm)...")
        hhs_start_shaker(ham_int, HHS3_MIDI.node, 1000)
        shake_timer = start_timer(10)
        shake_timer.wait(skip=device_simulation)
        hhs_stop_shaker(ham_int, HHS3_MIDI.node)
        
        # Transport to magnet
        print("Step 5: Moving to magnet...")
        transport_resource(ham_int, HHS3_MIDI.resource, MIDI_OnMagnet,
                         resource_type=GrippedResource.MIDI, core_gripper=True)
        
        # Let beads settle
        print("Step 6: Settling beads (60 sec)...")
        settle_timer = start_timer(60)
        settle_timer.wait(skip=device_simulation)
        
        # Remove supernatant
        print("Step 7: Removing supernatant...")
        double_aspirate_supernatant_96(ham_int, tracked_tips_300uL, tip_support, num_samples,
                                     source_plate=MIDI_OnMagnet, destination_plate=MPH_Waste,
                                     first_volume=200, second_volume=50, liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                                     aspiration_height=0, dispense_height=5)
        
        # Ethanol wash (2x) - Round 1
        print("Step 8: Performing ethanol washes (Round 1)...")
        for wash_num in range(2):
            print(f"  Ethanol wash {wash_num + 1}/2...")
            ethanol_wash(ham_int, tracked_tips_300uL, tip_support, num_samples,
                        ethanol_plate=EthanolReservoir, magnet_plate=MIDI_OnMagnet,
                        waste_plate=MIDI_Waste, wash_volume=post_shear_etoh_wash_volume,
                        first_removal_volume=200, second_removal_volume=50,
                        liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')
        
        # Air dry
        print("Step 9: Air drying beads (60 sec)...")
        dry_timer = start_timer(60)
        dry_timer.wait(skip=device_simulation)
        
        # Add Nuclease-Free Water
        print("Step 10: Adding elution buffer...")
        MIDI_OnMagnet_positions = [(MIDI_OnMagnet, idx) for idx in range(num_samples)]
        pip_transfer(ham_int, tracked_tips_300uL, Nuclease_Free_Water_positions,
                    MIDI_OnMagnet_positions, [post_shear_elution_buffer_volume]*num_samples,
                    liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
                    aspiration_height=1, dispense_height=1)
        
        # Transport to HHS3 to shake and incubate
        print("Step 11: Mixing elution on HHS3...")
        transport_resource(ham_int, MIDI_OnMagnet, HHS3_MIDI.resource,
                         resource_type=GrippedResource.MIDI, core_gripper=True)
        
        hhs_start_shaker(ham_int, HHS3_MIDI.node, 1000)
        shake_timer = start_timer(30)
        shake_timer.wait(skip=device_simulation)
        hhs_stop_shaker(ham_int, HHS3_MIDI.node)
        
        # Transport to magnet
        print("Step 12: Moving to magnet...")
        transport_resource(ham_int, HHS3_MIDI.resource, MIDI_OnMagnet,
                         resource_type=GrippedResource.MIDI, core_gripper=True)
        
        # Let beads settle
        print("Step 13: Settling beads (60 sec)...")
        settle_timer = start_timer(60)
        settle_timer.wait(skip=device_simulation)

        # Get fresh HSP plate from stack
        print("Step 14: Getting fresh HSP plate from stack...")
        transport_resource(ham_int, HSP_Stack.fetch_next(), HSP_Pipette2,
                             resource_type=GrippedResource.PCR, core_gripper=True)


        # Remove eluted samples
        print("Step 15: Removing eluted samples...")
        double_aspirate_supernatant_96(ham_int, tracked_tips_300uL, tip_support, num_samples,
                                     source_plate=MIDI_OnMagnet, destination_plate=HSP_Pipette2,
                                     first_volume=post_shear_elution_volume - 5, second_volume=5,
                                     liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                                     second_aspiration_height=0.7, dispense_height=1)

if __name__ == "__main__":
    # Run with simulation mode
    # Set simulating=False to run on real hardware
    # Set device_simulation=False to use real HHS devices
    magnetic_bead_cleanup(simulating=True, device_simulation=True)