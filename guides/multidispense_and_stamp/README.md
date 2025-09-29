# Multi-Dispense and Stamp

NGS protocols often require preparing master mixes by combining multiple reagents, then distributing them across samples. This example shows how to use multi-dispense for efficient reagent distribution and 96-channel transfer for high-throughput stamping.

## Setup Resources and Reagents

First, set up your deck layout with reagent carriers, plates, and tip racks:

```python
lmgr = LayoutManager(get_parent_lay_file())

CAR_VIALS_SMALL = layout_item(lmgr, ReagentTrackedEppiCarrier32, 'CAR_VIALS_SMALL')
CPAC_Reagents = layout_item(lmgr, ReagentTrackedPlate96, 'CPAC_HSP_0001')
MIDI_CPAC = layout_item(lmgr, Plate96, 'MIDI_CPAC')
MIDI_Plate = layout_item(lmgr, Plate96, 'MIDI_Plate')
HSP_Plate2 = layout_item(lmgr, Plate96, 'HSP_Plate2')
```

Assign reagents to their positions on carriers:

```python
fragmentation_buffer_positions = CAR_VIALS_SMALL.assign_reagent_map('FragmentationBuffer', [0])
buffer_eb_positions = CAR_VIALS_SMALL.assign_reagent_map('BufferEB', [2])
fragmentation_enzyme_positions = CPAC_Reagents.assign_reagent_map('FragmentationEnzyme', [0])
```

## Mix Source Reagents

Before multi-dispensing, mix reagents that may have settled. This is especially important for viscous reagents or those stored at room temperature:

```python
pip_mix(ham_int, tracked_tips_50uL, fragmentation_buffer_positions,
        mix_volume=200, mix_cycles=10,
        liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty')
```

## Multi-Dispense Master Mix Components

Multi-dispense efficiently distributes reagents from a single source to multiple destinations. This is ideal for adding the same reagent to many samples.

First, define your destination positions:

```python
num_samples = 32
mix_position = [(MIDI_CPAC, idx) for idx in range(num_samples)]
```

Calculate volumes based on your sample count:

```python
buffer_eb_vol = 20 * num_samples
frag_buffer_vol = 5 * num_samples
frag_enzyme_vol = 10 * num_samples
```

Add each reagent using multi-dispense. This aspirates once from the source and dispenses to all destinations:

```python
# Add Buffer EB (room temp)
multi_dispense(ham_int, tips=tracked_tips_300uL, source_positions=buffer_eb_positions,
               dispense_positions=mix_position, volumes=[buffer_eb_vol] * num_samples,
               aspiration_height=0, liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')

# Add Fragmentation Buffer (room temp)
multi_dispense(ham_int, tips=tracked_tips_300uL, source_positions=fragmentation_buffer_positions,
               dispense_positions=mix_position, volumes=[frag_buffer_vol] * num_samples,
               aspiration_height=0, liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')

# Add Fragmentation Enzyme (cold reagent from CPAC)
multi_dispense(ham_int, tracked_tips_300uL, fragmentation_enzyme_positions, mix_position,
               volumes=[frag_enzyme_vol] * num_samples, aspiration_height=0,
               liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
```

## Mix the Master Mix

After adding all components, mix the combined reagents:

```python
pip_mix(ham_int, tips=tracked_tips_300uL, positions_to_mix=mix_position,
        mix_cycles=3, mix_volume=25, liquid_height=0,
        liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')
```

## Transport Plate

Move the MIDI plate from the CPAC (temperature-controlled carrier) to the deck for stamping. Use the core gripper for MIDI plate transport:

```python
transport_resource(ham_int, MIDI_CPAC, MIDI_Plate,
                   resource_type=GrippedResource.MIDI, core_gripper=True)
```

## 96-Channel Transfer (Stamping)

Use the 96-channel head to transfer the master mix to your sample plate. This performs parallel transfers across all wells:

```python
transfer_96(ham_int, tracked_tips_300uL, tip_support, num_samples,
            source_plate=MIDI_Plate, target_plate=HSP_Plate2, volume=40,
            liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty')
```

The `transfer_96` function handles tip pickup from the support, aspiration from the source plate, dispensing to the target plate, and tip ejection.