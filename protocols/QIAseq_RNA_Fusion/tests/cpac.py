#!/usr/bin/env python3
"""CPAC Temperature Control Test"""

from pyhamilton import HamiltonInterface, LayoutManager, layout_item
from pyhamilton.devices import (initialize_cpac, set_temperature_target_cpac, 
                                start_temperature_control_cpac, get_temperature_cpac)
from pyhamilton.resources import Plate96
from pyhamilton.transport import transport_resource, GrippedResource
import time

def test_cpac():
    lmgr = LayoutManager('QIAseq RNA Fusion XP Panels.lay')
    HSP_CPAC = layout_item(lmgr, Plate96, 'HSP_CPAC')
    HSP_Pipette = layout_item(lmgr, Plate96, 'HSP_Pipette')
    
    with HamiltonInterface(simulating=True, windowed=True) as ham_int:
        ham_int.initialize()
        
        # Initialize and set temperature
        initialize_cpac(ham_int, controller_id=1, simulating=True)
        set_temperature_target_cpac(ham_int, controller_id=1, device_id=1, target_temp=4.0)
        start_temperature_control_cpac(ham_int, controller_id=1, device_id=1)
        
        # Test temperature reading
        temp = get_temperature_cpac(ham_int, controller_id=1, device_id=1)
        print(f"CPAC Temperature: {temp}°C")
        
        # Test plate transport
        transport_resource(ham_int, HSP_Pipette, HSP_CPAC, 
                            resource_type=GrippedResource.PCR, core_gripper=True)
        transport_resource(ham_int, HSP_CPAC, HSP_Pipette,
                            resource_type=GrippedResource.PCR, core_gripper=True)
        
        print("CPAC: OK")
        

if __name__ == "__main__":
    test_cpac()