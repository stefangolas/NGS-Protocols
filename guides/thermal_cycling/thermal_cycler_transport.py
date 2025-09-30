"""Thermal Cycling Test"""

from pyhamilton import HamiltonInterface, LayoutManager, layout_item
from pyhamilton.ngs import protocol
from pyhamilton.resources import Plate96, Lid, StackedResources
from pyhamilton.transport import transport_resource, GrippedResource, GripDirection
from pyhamilton.devices import (odtc_connect, odtc_initialize, odtc_open_door, 
                                odtc_close_door, odtc_execute_protocol, odtc_wait_for_idle)
from pathlib import Path

def get_parent_lay_file():
    parent = Path(__file__).resolve().parent.parent  # go one level up
    for file in parent.glob("*.lay"):
        return str(file)  # return first match
    return None  # if no .lay file found

def thermal_cycling_with_transport(simulating):
    lmgr = LayoutManager(get_parent_lay_file())

    HSP_Pipette2 = layout_item(lmgr, Plate96, 'HSP_Pipette2')
    HSP_CPAC = layout_item(lmgr, Plate96, 'HSP_CPAC')
    HSP_ODTC = layout_item(lmgr, Plate96, 'HSP_ODTC')
    HSP_ODTC_Lid = layout_item(lmgr, Lid, 'Ham_ComfortLid_ODTC')
    
    Lid_Stack = StackedResources.from_prefix(
        tracker_id="Ham_ComfortLid_ParkPos", prefix="Ham_ComfortLid_ParkPos",
        count=4, lmgr=lmgr, resource_type=Lid
    )

    windowed = not simulating
    with HamiltonInterface(simulating=simulating, windowed=windowed) as ham_int:
        ham_int.initialize()
        
        # Initialize ODTC
        device_id = odtc_connect(ham_int, simulation_mode=True, 
                                local_ip='1.2.3.4', device_ip='5.6.7.8', device_port='COM4')
        odtc_initialize(ham_int, device_id=device_id)
                
        odtc_open_door(ham_int, device_id=device_id)

        # Move PCR plate to thermal cycler
        transport_resource(ham_int, HSP_Pipette2, HSP_ODTC,
                            grip_direction=GripDirection.RIGHT,
                            resource_type=GrippedResource.PCR, iswap=True)
        
        # Add lid to PCR plate from lid stack
        transport_resource(ham_int, Lid_Stack.fetch_next(), HSP_ODTC_Lid,
                            grip_direction=GripDirection.RIGHT,
                            resource_type=GrippedResource.LID, iswap=True)
        
        odtc_close_door(ham_int, device_id=device_id)
        
        # Run protocol
        odtc_execute_protocol(ham_int, device_id=device_id,
                            method_name='FirstStrandDNASynthesis.xml', simulating=True)
        
        # Wait for thermal cycler
        odtc_wait_for_idle(ham_int, device_id=device_id, simulating=simulating, check_interval=5)

        # Unload plate
        odtc_open_door(ham_int, device_id=device_id)
        transport_resource(ham_int, HSP_ODTC_Lid, Lid_Stack.put_back(),
                            grip_direction=GripDirection.RIGHT,
                            resource_type=GrippedResource.LID, iswap=True)
        transport_resource(ham_int, HSP_ODTC, HSP_CPAC,
                            grip_direction=GripDirection.RIGHT,
                            resource_type=GrippedResource.PCR, iswap=True)
        odtc_close_door(ham_int, device_id=device_id)
            

if __name__ == "__main__":
    thermal_cycling_with_transport(simulating=False)