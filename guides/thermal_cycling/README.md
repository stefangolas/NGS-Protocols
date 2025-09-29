# Thermal Cycling with Transport

NGS protocols often involve using the thermal cycler to perform reactions like end-tailing and adapter ligation.

The ODTC should be connected to and initialized at the beginning of the protocol.
```python
device_id = odtc_connect(ham_int, simulation_mode=True, 
                        local_ip='1.2.3.4', device_ip='5.6.7.8', device_port='COM4')
odtc_initialize(ham_int, device_id=device_id)
        
```

Then we move the desired PCR plate from the deck to the ODTC using the iSWAP. This always uses a right-side grip
with PCR plate gripper settings.
```python
odtc_open_door(ham_int, device_id=device_id)

# Move PCR plate to thermal cycler
transport_resource(ham_int, HSP_Pipette2, HSP_ODTC,
                    grip_direction=GripDirection.RIGHT,
                    resource_type=GrippedResource.PCR, iswap=True)
```

Move a lid onto the PCR plate from the lid stack. This uses pre-defined gripper settings for lids.
```python
transport_resource(ham_int, Lid_Stack.fetch_next(), HSP_ODTC_Lid,
                    grip_direction=GripDirection.RIGHT,
                    resource_type=GrippedResource.LID, iswap=True)
```

Close the door, execute the protocol, and call the wait function to block the execution of the rest of the protocol.
```python
odtc_close_door(ham_int, device_id=device_id)

# Run protocol
odtc_execute_protocol(ham_int, device_id=device_id,
                    method_name='FirstStrandDNASynthesis.xml', simulating=True)

# Wait for thermal cycler
odtc_wait_for_idle(ham_int, device_id=device_id, simulating=simulating, check_interval=5)
```

Once the ODTC is finished, we open the door and move the PCR plate back to the deck.
```python
# Unload plate
odtc_open_door(ham_int, device_id=device_id)
transport_resource(ham_int, HSP_ODTC_Lid, Lid_Stack.put_back(),
                    grip_direction=GripDirection.RIGHT,
                    resource_type=GrippedResource.LID, iswap=True)
transport_resource(ham_int, HSP_ODTC, HSP_CPAC,
                    grip_direction=GripDirection.RIGHT,
                    resource_type=GrippedResource.PCR, iswap=True)
odtc_close_door(ham_int, device_id=device_id)
```