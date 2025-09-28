from pyhamilton import (HamiltonInterface, LayoutManager, Plate96, Tip96, normal_logging, layout_item)

from pyhamilton.resources import Lid, StackedResources
from pyhamilton.transport import (transport_resource, GripDirection, GrippedResource)
from pyhamilton.devices import odtc_open_door, odtc_close_door, odtc_execute_protocol, odtc_get_status, odtc_initialize, odtc_connect

import time
import os


def thermal_cycle_with_plate_movement(ham_int, source_plate_sequence, odtc_location_sequence, lid_stack, protocol_file_path):
    """
    Common pattern for thermal cycling: Move plate to thermal cycler, add lid, run program, remove lid, move back
    
    Args:
        ham_int: Hamilton interface object
        source_plate_name: Name of the source plate position (e.g., 'HSP_Pipette2')
        destination_plate_name: Name of the destination plate position (e.g., 'HSP_ODTC')
        lid_stack: StackedResources object for lids
        thermal_cycler: ODTC thermal cycler object
    """
    # Open thermal cycler door
    odtc_open_door(ham_int, device_id)

    # Move plate from source to destination (usually ODTC) using iswap
    transport_resource(ham_int, source_plate_sequence, odtc_location_sequence, 
                      grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
    
    # Add lid to the plate
    lid = lid_stack.fetch_next()
    transport_resource(ham_int, lid, 'Ham_ComfortLid_ODTC', 
                      resource_type=GrippedResource.LID, grip_direction=GripDirection.RIGHT, iswap=True)

    # Close door and run thermal cycler program
    odtc_close_door(ham_int, device_id)
    duration, resultID = odtc_execute_protocol(ham_int, device_id, protocol_file_path, priority=1)
    odtc_ready = False
    while not odtc_ready:
        odtc_ready = odtc_get_status(ham_int, device_id) == 'idle'
        time.sleep(5)

    
    odtc_get_status(ham_int, device_id)
    
    # Open door after program completion
    odtc_open_door(ham_int, device_id)

    # Remove lid from the plate
    transport_resource(ham_int, 'Ham_ComfortLid_ODTC', lid_stack.fetch_next_unoccupied(), 
                      resource_type=GrippedResource.LID, core_gripper=True)
    
    # Move plate back to source position using iswap
    transport_resource(ham_int, odtc_location_sequence, source_plate_sequence, 
                      resource_type=GrippedResource.PCR, iswap=True)
    
    return duration, resultID


lmgr = LayoutManager('LSK109_MPH_OxfordNanopore_v0.9.2.lay')

LidStack = StackedResources.from_prefix(
    tracker_id="Ham_ComfortLid_Stack_ParkPos",
    prefix="Ham_ComfortLid_Stack_ParkPos",
    count=3,
    lmgr=lmgr,
    resource_type=Lid,
    reset=True)

HSP_Pipette2 = layout_item(lmgr, Plate96, 'HSP_Pipette2')
HSP_ODTC = layout_item(lmgr, Plate96, 'HSP_ODTC')
odtc_lid = layout_item(lmgr, Lid, 'Ham_ComfortLid_ODTC')

if __name__ == "__main__":
    with HamiltonInterface(windowed=True, simulating=False) as ham_int:
        ham_int.initialize()
        normal_logging(ham_int, os.getcwd())

        # Create thermal cycler instance
        device_id = odtc_connect(ham_int, simulation_mode=True, local_ip='192.168.1.200', device_ip='192.168.1.50')
        odtc_initialize(ham_int, device_id)
        odtc_open_door(ham_int, device_id)
        transport_resource(ham_int, HSP_Pipette2, HSP_ODTC, grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
        
        # Move lid from lid stack to ODTC
        lid = LidStack.fetch_next()
        print("Lid fetched from stack:", lid)
        transport_resource(ham_int, lid, odtc_lid, 
                      resource_type=GrippedResource.LID, grip_direction=GripDirection.RIGHT, iswap=True)

        response = odtc_execute_protocol(ham_int, device_id, 'protocol.xml', simulating=True, priority=1)
        odtc_ready = False
        while not odtc_ready:
            response = odtc_get_status(ham_int, device_id, simulating=True)
            odtc_ready = response.state == 'idle'
            time.sleep(5)

        
        odtc_get_status(ham_int, device_id, simulating=True)
        
        # Open door after program completion
        odtc_open_door(ham_int, device_id)

        transport_resource(ham_int, odtc_lid, lid, 
                      resource_type=GrippedResource.LID, grip_direction=GripDirection.RIGHT, iswap=True)

        transport_resource(ham_int, HSP_ODTC, HSP_Pipette2, grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)

