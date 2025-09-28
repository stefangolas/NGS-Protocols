#!/usr/bin/env python3
"""HHS Device Test - Tests all 5 Hamilton Heater Shaker nodes"""

from pyhamilton import HamiltonInterface, hhs_set_simulation
from pyhamilton.devices import hhs_create_usb_device, hhs_start_shaker, hhs_stop_shaker
import time

    

def test_hhs_devices():
    with HamiltonInterface(simulating=True, windowed=True) as ham_int:
        ham_int.initialize()
        hhs_set_simulation(ham_int, 1)
        
        # Test all 5 HHS nodes
        for node in range(1, 6):
            try:
                hhs_create_usb_device(ham_int, node)
                hhs_start_shaker(ham_int, node, speed=1000)
                time.sleep(1)
                hhs_stop_shaker(ham_int, node)
                print(f"HHS {node}: OK")
            except Exception as e:
                print(f"HHS {node}: FAILED - {e}")

if __name__ == "__main__":
    test_hhs_devices()