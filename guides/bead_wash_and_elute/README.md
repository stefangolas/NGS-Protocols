# Magnetic Bead Wash and Elute

Many NGS protocols involve using magnetic beads to wash sample mixtures, especially
after performing a reaction step in a thermal cycler. We add a magnetic bead
solution to the sample mixture, the DNA or RNA product binds to the beads during an
incubation, and then we move the sample plate to a magnet stand. We can aspirate supernatant
with dissolved impurities while the beads are pulled down by the magnet, then move the plate
off the stand to resuspend beads in wash buffer, and repeat until the desired purity
is acheived.

This process can also be used to perform DNA size selection as explained here.

Note that supernatant aspirations are performed with the `double_aspirate_supernatant_96(...)`
function. This performs two sequential aspirations of supernatant. The first uses capacitive liquid level
detection and liquid following to aspirate most of the supernatant, 
and the second aspirates a typically smaller volume from a fixed height. This allows residual 
supernatant to settle from the walls of the well and ensures that the bead is not disturbed on the second aspiration.  

Get a MIDI plate from the stack, and add magnetic beads and reaction mixture to it.
```python
# Get fresh MIDI plate from stack and move to HHS3
print("Step 1: Preparing MIDI plate on HHS3...")
transport_resource(ham_int, MIDI_Stack.fetch_next(), HHS3_MIDI.resource,
                    resource_type=GrippedResource.MIDI, core_gripper=True)

# Add QIAseq Beads to HHS3 position (Round 1)
print("Step 2: Adding QIAseq Beads (Round 1)...")
HHS3_MIDI_positions = [(HHS3_MIDI.resource, idx) for idx in range(num_samples)]
volumes = [post_shear_magbead_volume] * num_samples
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
```

Move MIDI plate to magnet, wait for beads to settle, aspirate off reaction mixture supernatant, then perform ethanol washes.
```python
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
```

Dry beads (be careful not to overdry), then add nuclease-free water and incubate. This elutes the DNA or RNA product off of the beads
and into solution. Carefully transfer supernatant to fresh HSP plate.
```python
# Air dry
print("Step 9: Air drying beads (60 sec)...")
dry_timer = start_timer(60)
dry_timer.wait(skip=device_simulation)

# Add Nuclease-Free Water
print("Step 10: Adding elution buffer (Round 1)...")
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
```