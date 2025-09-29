# Magnetic Bead Wash and Elute

Many NGS protocols involve using magnetic beads to purify DNA or RNA products after reactions. The beads bind nucleic acids during incubation, allowing you to wash away impurities while the beads are held by a magnet, then elute the purified product in clean buffer.

This tutorial covers the complete bead cleanup workflow: adding beads, binding samples, washing with ethanol, and eluting purified product.

## Overview of the Process

1. Prepare a fresh MIDI plate and add magnetic beads
2. Transfer reaction mixture to the beads and mix
3. Move plate to magnet and remove supernatant
4. Perform ethanol washes while beads are held by magnet
5. Air dry the beads
6. Add elution buffer and incubate
7. Collect purified eluate

## Step 1: Get Fresh MIDI Plate

Start by getting a fresh MIDI plate from the stack and placing it on the heated shaker (HHS3):

```python
print("Step 1: Preparing MIDI plate on HHS3...")
transport_resource(ham_int, MIDI_Stack.fetch_next(), HHS3_MIDI.resource,
                   resource_type=GrippedResource.MIDI, core_gripper=True)
```

The `MIDI_Stack.fetch_next()` method automatically tracks which plate to pick next from the stack.

## Step 2: Add Magnetic Beads

Add the magnetic bead reagent to each sample position on the MIDI plate:

```python
print("Step 2: Adding QIAseq Beads (Round 1)...")
# Define destination positions for all samples
HHS3_MIDI_positions = [(HHS3_MIDI.resource, idx) for idx in range(num_samples)]

# Calculate volumes for each sample
volumes = [post_shear_magbead_volume] * num_samples

# Transfer beads to each position
pip_transfer(ham_int, tracked_tips_300uL, QIAseq_Beads_positions,
             HHS3_MIDI_positions, volumes,
             liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
             aspiration_height=1, dispense_height=1)
```

The bead volume is typically calculated based on the sample volume and desired binding ratio (e.g., 1.8x sample volume for size selection).

## Step 3: Transfer Reaction Samples

Use the 96-channel head to transfer your reaction samples from the thermal cycler output plate to the MIDI plate containing beads:

```python
print("Step 3: Transferring samples to MIDI plate...")
transfer_96(ham_int, tracked_tips_300uL, tip_support=tip_support, num_samples=num_samples,
            target_plate=HSP_CPAC, source_plate=HHS3_MIDI.resource, volume=sample_volume,
            aspiration_mix_cycles=3, aspiration_mix_volume=sample_volume,
            liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
            aspiration_height=1, dispense_height=1)
```

Note the `aspiration_mix_cycles=3` parameter - this mixes the samples before aspirating to ensure homogeneous transfer.

## Step 4: Mix Beads with Samples

Mix the beads thoroughly with the samples using the heated shaker:

```python
print("Step 4: Mixing on shaker (10 sec @ 1000 rpm)...")
# Start the shaker at 1000 rpm
hhs_start_shaker(ham_int, HHS3_MIDI.node, 1000)

# Create a timer for the mixing duration
shake_timer = start_timer(10)
shake_timer.wait(skip=device_simulation)

# Stop the shaker
hhs_stop_shaker(ham_int, HHS3_MIDI.node)
```

This rapid mixing ensures the nucleic acids bind efficiently to the magnetic beads.

## Step 5: Move Plate to Magnet

Transport the MIDI plate from the shaker to the magnetic plate stand:

```python
print("Step 5: Moving to magnet...")
transport_resource(ham_int, HHS3_MIDI.resource, MIDI_OnMagnet,
                   resource_type=GrippedResource.MIDI, core_gripper=True)
```

## Step 6: Let Beads Settle

Allow time for the magnetic beads to pellet against the side of the wells:

```python
print("Step 6: Settling beads (60 sec)...")
settle_timer = start_timer(60)
settle_timer.wait(skip=device_simulation)
```

The settling time depends on bead type and volume - typically 30-120 seconds.

## Step 7: Remove Reaction Supernatant

Use the double aspirate function to remove the reaction supernatant without disturbing the bead pellet:

```python
print("Step 7: Removing supernatant...")
double_aspirate_supernatant_96(ham_int, tracked_tips_300uL, tip_support, num_samples,
                               source_plate=MIDI_OnMagnet, destination_plate=MPH_Waste,
                               first_volume=200, second_volume=50, 
                               liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                               aspiration_height=0, dispense_height=5)
```

The `double_aspirate_supernatant_96()` function performs two sequential aspirations:
- **First aspiration**: Uses capacitive liquid level detection (CLLD) to find the liquid surface and aspirates most of the supernatant while following the liquid level down
- **Second aspiration**: Removes residual supernatant from a fixed height after it has settled from the well walls

This two-step approach maximizes supernatant removal while minimizing the risk of disturbing the magnetic bead pellet.

## Step 8: Perform Ethanol Washes

Wash the beads twice with ethanol to remove contaminants:

```python
print("Step 8: Performing ethanol washes (Round 1)...")
for wash_num in range(2):
    print(f"  Ethanol wash {wash_num + 1}/2...")
    ethanol_wash(ham_int, tracked_tips_300uL, tip_support, num_samples,
                 ethanol_plate=EthanolReservoir, magnet_plate=MIDI_OnMagnet,
                 waste_plate=MIDI_Waste, wash_volume=post_shear_etoh_wash_volume,
                 first_removal_volume=200, second_removal_volume=50,
                 liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')
```

Each `ethanol_wash()` call:
1. Adds ethanol to the beads
2. Waits briefly for mixing by diffusion
3. Removes the ethanol using double aspiration

The beads remain on the magnet throughout the wash steps.

## Step 9: Air Dry Beads

Dry the beads to remove residual ethanol:

```python
print("Step 9: Air drying beads (60 sec)...")
dry_timer = start_timer(60)
dry_timer.wait(skip=device_simulation)
```

**Important**: Be careful not to over-dry the beads. Over-dried beads are difficult to resuspend and can reduce yield. Beads should appear matte but not cracked.

## Step 10: Add Elution Buffer

Add nuclease-free water or elution buffer to the dried beads:

```python
print("Step 10: Adding elution buffer (Round 1)...")
# Define positions while plate is still on magnet
MIDI_OnMagnet_positions = [(MIDI_OnMagnet, idx) for idx in range(num_samples)]

# Add elution buffer to beads
pip_transfer(ham_int, tracked_tips_300uL, Nuclease_Free_Water_positions,
             MIDI_OnMagnet_positions, [post_shear_elution_buffer_volume]*num_samples,
             liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty',
             aspiration_height=1, dispense_height=1)
```

The elution volume is typically much smaller than the original sample volume to concentrate the product.

## Step 11: Mix and Incubate for Elution

Move the plate off the magnet to resuspend the beads and elute the nucleic acids:

```python
print("Step 11: Mixing elution on HHS3...")
# Move from magnet to shaker
transport_resource(ham_int, MIDI_OnMagnet, HHS3_MIDI.resource,
                   resource_type=GrippedResource.MIDI, core_gripper=True)

# Mix to resuspend beads and elute product
hhs_start_shaker(ham_int, HHS3_MIDI.node, 1000)
shake_timer = start_timer(30)
shake_timer.wait(skip=device_simulation)
hhs_stop_shaker(ham_int, HHS3_MIDI.node)
```

This mixing step is critical - the nucleic acids need to be released from the beads into the elution buffer.

## Step 12: Return to Magnet

Move the plate back to the magnet to pellet the beads again:

```python
print("Step 12: Moving to magnet...")
transport_resource(ham_int, HHS3_MIDI.resource, MIDI_OnMagnet,
                   resource_type=GrippedResource.MIDI, core_gripper=True)
```

## Step 13: Let Beads Settle Again

Allow the beads to pellet so you can cleanly transfer the eluate:

```python
print("Step 13: Settling beads (60 sec)...")
settle_timer = start_timer(60)
settle_timer.wait(skip=device_simulation)
```

## Step 14: Prepare Destination Plate

Get a fresh HSP plate from the stack to receive the purified eluate:

```python
print("Step 14: Getting fresh HSP plate from stack...")
transport_resource(ham_int, HSP_Stack.fetch_next(), HSP_Pipette2,
                   resource_type=GrippedResource.PCR, core_gripper=True)
```

## Step 15: Transfer Purified Eluate

Carefully transfer the eluate (containing purified nucleic acids) to the fresh plate:

```python
print("Step 15: Removing eluted samples...")
double_aspirate_supernatant_96(ham_int, tracked_tips_300uL, tip_support, num_samples,
                               source_plate=MIDI_OnMagnet, destination_plate=HSP_Pipette2,
                               first_volume=post_shear_elution_volume - 5, second_volume=5,
                               liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty',
                               second_aspiration_height=0.7, dispense_height=1)
```

Note the careful volume splitting:
- First aspiration: Takes most of the elution volume
- Second aspiration: Takes the final 5 µL from a slightly elevated height (0.7 mm) to avoid aspirating any beads

The `second_aspiration_height=0.7` parameter is critical here - it ensures you stay above the bead pellet while maximizing eluate recovery.

## Key Parameters to Optimize

When adapting this protocol, consider adjusting:
- **Bead volume**: Ratio to sample volume affects binding capacity and size selection
- **Settling times**: Longer for larger volumes or viscous samples
- **Elution volume**: Smaller volumes give higher concentration but lower total yield
- **Elution time**: Longer incubation can improve yield
- **Drying time**: Must be sufficient to remove ethanol but not crack the beads