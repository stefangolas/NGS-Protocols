# Deck Setup and Resource Management

Before running any Hamilton liquid handling protocol, you need to set up your deck layout with all the plates, tip racks, reagent containers, and other resources. This tutorial shows how to properly initialize and organize your deck for a complex NGS library prep workflow.

## Overview

A complete deck setup involves:
1. Loading the deck layout file
2. Defining all labware positions
3. Mapping reagents to their containers
4. Setting up stacked resources (plate stacks, tip racks)
5. Configuring tracked tips for liquid classes
6. Initializing the Hamilton interface

## Step 1: Load the Deck Layout File

Every Hamilton protocol needs a `.lay` layout file that defines the physical positions of all labware on the deck:

```python
def get_parent_lay_file():
    """Find the .lay file in the parent directory."""
    parent = Path(__file__).resolve().parent.parent
    for file in parent.glob("*.lay"):
        return str(file)
    return None
```

Load this layout file with the LayoutManager:

```python
lay_file = get_parent_lay_file()
lmgr = LayoutManager(lay_file)
```

The `LayoutManager` reads the layout file and provides access to all defined positions by name.

## Step 2: Define Protocol Parameters

Set key parameters at the beginning to make your protocol easily adjustable:

```python
num_samples = 8
```

This parameter controls how many samples will be processed and is used throughout the protocol.

## Step 3: Map PCR Plate Positions

Define all your PCR plate positions for thermal cycling and standard pipetting operations:

```python
# PCR plates for thermal cycling steps
HSP_Plate = layout_item(lmgr, Plate96, 'HSP_Pipette')
HSP_Plate2 = layout_item(lmgr, Plate96, 'HSP_Pipette2')
HSP_Waste = layout_item(lmgr, Plate96, 'HSP_Waste')
HHS1_HSP = layout_item(lmgr, Plate96, 'HHS1_HSP')
HSP_Magnet = layout_item(lmgr, Plate96, 'HSP_OnMagnet')
HSP_ODTC = layout_item(lmgr, Plate96, 'HSP_ODTC')
HSP_ODTC_Lid = layout_item(lmgr, Lid, 'Ham_ComfortLid_ODTC')
```

The `layout_item()` function maps a physical position name from the `.lay` file to a resource type:
- First argument: The layout manager
- Second argument: The resource class (e.g., `Plate96`, `Lid`)
- Third argument: The position name exactly as defined in the layout file

**Common PCR plate positions:**
- **HSP_Plate / HSP_Plate2**: General pipetting positions
- **HSP_Waste**: Waste collection
- **HHS1_HSP**: Heated shaker position for incubations
- **HSP_Magnet**: Magnetic plate stand for bead separations
- **HSP_ODTC**: Inside the thermal cycler
- **HSP_ODTC_Lid**: Lid staging position for thermal cycler

## Step 4: Map MIDI Plate Positions

MIDI plates are deeper than standard PCR plates, making them ideal for magnetic bead work where higher volumes are needed:

```python
# MIDI plates for magnetic bead operations
MIDI_Plate = layout_item(lmgr, Plate96, 'MIDI_Pipette')
MIDI_Waste = layout_item(lmgr, Plate96, 'MIDI_Waste')
MIDI_CPAC = layout_item(lmgr, Plate96, 'MIDI_CPAC')
MIDI_OnMagnet = layout_item(lmgr, Plate96, 'MIDI_OnMagnet')
HHS5_MIDI = layout_item(lmgr, Plate96, 'HHS5_MIDI')
```

**Common MIDI plate positions:**
- **MIDI_Plate**: General pipetting position
- **MIDI_CPAC**: Temperature-controlled carrier (cooling) position
- **MIDI_OnMagnet**: MIDI-sized magnetic plate stand
- **HHS5_MIDI**: Heated shaker for MIDI plates

## Step 5: Define Reagent Containers

Map all reagent containers to their deck positions. Different container types are used based on volume and temperature requirements:

```python
# Small volume reagents (tubes in carrier)
CAR_VIALS_SMALL = layout_item(lmgr, ReagentTrackedEppiCarrier32, 'CAR_VIALS_SMALL')

# Cold-sensitive reagents (plate on CPAC)
CPAC_Reagents = layout_item(lmgr, ReagentTrackedPlate96, 'CPAC_HSP_0001')

# Medium bulk reagents (60 mL reservoirs)
RGT_01 = layout_item(lmgr, ReagentTrackedReservoir60mL, 'rgt_cont_60ml_BC_A00_0001')
RGT_02 = layout_item(lmgr, ReagentTrackedReservoir60mL, 'rgt_cont_60ml_BC_A00_0002')

# Large bulk reagents (bulk plate for ethanol washes)
EthanolReservoir = layout_item(lmgr, ReagentTrackedBulkPlate, 'RGT_Ethanol')
```

**Choosing the right container type:**
- `ReagentTrackedEppiCarrier32`: For small volumes (< 2 mL) in tubes
- `ReagentTrackedPlate96`: For reagents stored in wells on a plate, like on a CPAC adapter.
- `ReagentTrackedReservoir60mL`: For medium bulk reagents (up to 60 mL)
- `ReagentTrackedBulkPlate`: For high-volume bulk reagents like ethanol

All `ReagentTracked*` classes automatically track liquid levels and volumes dispensed.

## Step 6: Assign Reagents to Positions

Use the `assign_reagent_map()` method to tell the system which reagent is in which position:

```python
# Room temperature small volume reagents
fragmentation_buffer_positions = CAR_VIALS_SMALL.assign_reagent_map('FragmentationBuffer', [0])
ligation_mix_positions = CAR_VIALS_SMALL.assign_reagent_map('LigationMix', [1])
buffer_eb_positions = CAR_VIALS_SMALL.assign_reagent_map('BufferEB', [2])
```

The `assign_reagent_map()` method takes:
- **Reagent name**: A string identifier for the reagent
- **Position indices**: A list of positions (0-indexed) where this reagent is located

For cold-sensitive reagents on the CPAC:

```python
# Cold-sensitive reagents in CPAC
fragmentation_enzyme_positions = CPAC_Reagents.assign_reagent_map('FragmentationEnzyme', [0])
dna_ligase_positions = CPAC_Reagents.assign_reagent_map('DNALigase', [1])
library_amp_mix_positions = CPAC_Reagents.assign_reagent_map('LibraryAmpMix', [2])
```

For bulk reagents spanning multiple wells:

```python
# Bulk reagents spanning multiple reservoir positions
nuclease_free_water_positions = RGT_01.assign_reagent_map('NucleaseFreeWater', range(8))
spriselect_positions = RGT_02.assign_reagent_map('SPRIselect', range(8))
```

`range(8)` is always used with 60mL troughs because it's a single contiguous container.

## Step 7: Define Sample Positions

For reagents that are pre-loaded per sample (like indexing primers), define positions directly:

```python
# Index positions (one per sample)
index_positions = [(HHS1_HSP, i) for i in range(num_samples)]
```

This creates a list of tuples: `[(plate, well_index), ...]` for each sample. These are ready to be passed as a position variable to a liquid-handling command.

## Step 8: Set Up Liquid Waste

Define a position for liquid waste disposal:

```python
Liquid_Waste = layout_item(lmgr, Waste96, 'core96externalwaste_0001')
```

## Step 9: Configure Stacked Resources

Stacked resources automatically track which item to pick next from a stack of plates or lids:

```python
# HSP plate stack (5 plates)
HSP_Stack = StackedResources.from_prefix(
    tracker_id="BioRadHardShell_Stack4",
    prefix="BioRadHardShell_Stack4",
    lmgr=lmgr,
    resource_type=Plate96,
    count=5)
```

The `from_prefix()` method requires:
- **tracker_id**: Unique identifier for tracking
- **prefix**: Base name of positions in the layout file (e.g., "BioRadHardShell_Stack4_1", "BioRadHardShell_Stack4_2", etc.)
- **lmgr**: Layout manager instance
- **resource_type**: Type of resource in the stack
- **count**: Number of items in the stack

Set up lid stacks:

```python
# Lid stack (3 lids)
Lid_Stack = StackedResources.from_prefix(
    tracker_id="Ham_ComfortLid_Stack",
    prefix="Ham_ComfortLid_Stack",
    lmgr=lmgr,
    resource_type=Lid,
    count=3)
```

Set up MIDI plate stacks:

```python
# MIDI plate stack (3 plates)
MIDI_Stack = StackedResources.from_prefix(
    tracker_id="AbgeneMIDI_Stack1",
    prefix="AbgeneMIDI_Stack1",
    lmgr=lmgr,
    resource_type=Plate96,
    count=3)
```

**Using stacked resources in protocols:**
- `Stack.fetch_next()`: Gets the next available item from the stack
- `Stack.put_back()`: Returns an item to the stack

## Step 10: Set Up Tracked Tips

Tracked tips monitor tip usage and automatically advance to fresh tips:

```python
# 50 µL filter tips (8 racks)
tracked_tips_50uL = TrackedTips.from_prefix(
    tracker_id="TIP_50uLF_L",
    volume_capacity=50,
    prefix="TIP_50uLF_L",
    count=8,
    tip_type=Tip96,
    lmgr=lmgr)
```

Parameters:
- **tracker_id**: Unique identifier
- **volume_capacity**: Maximum volume for these tips (µL)
- **prefix**: Base name in layout file
- **count**: Number of tip racks
- **tip_type**: Resource type (usually `Tip96`)
- **lmgr**: Layout manager

Set up multiple tip sizes for different volumes:

```python
# 300 µL tips (8 racks) - most common for NGS
tracked_tips_300uL = TrackedTips.from_prefix(
    tracker_id="STF_L",
    volume_capacity=300,
    prefix="STF_L",
    count=8,
    tip_type=Tip96,
    lmgr=lmgr)

# 1000 µL tips (2 racks) - for bulk transfers
tracked_tips_1000uL = TrackedTips.from_prefix(
    tracker_id="HTF_L",
    volume_capacity=1000,
    prefix="HTF_L",
    count=2,
    tip_type=Tip96,
    lmgr=lmgr)
```

## Step 11: Configure Tip Support

The tip support is used for 96-channel operations where we pick up partial racks of tips using an offset to conserve tip usage.

```python
tip_support_resource = layout_item(lmgr, Tip96, 'TipSupport_0001')
tip_support = TipSupportTracker(tip_support_resource)
```

The `TipSupportTracker` wraps the tip support position and manages tip presence tracking.

## Step 12: Initialize the Hamilton Interface

Finally, initialize the Hamilton system and run your protocol:

```python
windowed = not simulating
with HamiltonInterface(simulating=simulating, windowed=windowed) as ham_int:
    ham_int.initialize()
    normal_logging(ham_int, os.getcwd())
    
    print("Deck layout initialized successfully")
    
    # Your protocol steps go here...
```

**Key initialization parameters:**
- **simulating**: If `True`, runs in simulation without hardware
- **windowed**: If `True`, shows the Hamilton Run Control window. From here you can run in Venus simulation or live.

The `ham_int.initialize()` call:
- Initializes all connected hardware
- Homes the robotic arms
- Verifies deck layout matches the physical setup

The `normal_logging()` function enables detailed logging to help with troubleshooting.


## Best Practices

**Organization:**
- Group related positions together (PCR plates, MIDI plates, reagents, etc.)
- Use clear, descriptive variable names
- Keep reagent assignments near their container definitions

**Reagent tracking:**
- Always use `ReagentTracked*` classes for containers with reagents so you can calculate required volumes from simulation mode.


**Tip management:**
- Run in simulation mode to calculate tip usage before running a protocol

**Simulation:**
- Always test in simulation mode first (`simulating=True`)
- Verify all positions are correctly mapped before running on hardware
- Use `device_simulation=True` for HHS devices during development

**Error prevention:**
- Double-check position names match your `.lay` file exactly
- Verify resource types match the physical labware
- Ensure stack counts match the actual number of items loaded