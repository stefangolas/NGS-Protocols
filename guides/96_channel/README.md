# 96-Channel Tip Management with Automatic Volume Switching

When using the Hamilton 96-channel head, tips must be managed differently than with single-channel pipetting. The TipSupportTracker and TrackedTips system provides automatic tip management, including intelligent switching between different tip volumes.

## The Challenge

96-channel operations often require different tip sizes for different volumes:
- Small volumes (< 50 µL) → 50 µL tips
- Medium/large volumes (50-300 µL) → 300 µL tips

Manually swapping tip racks between operations would be tedious and error-prone. The tip support system handles this automatically.

## How It Works

The **TipSupportTracker** maintains a single "staging area" where one tip rack is loaded at a time. When you request tips of a different volume than what's currently loaded, the system:

1. Detects the volume mismatch
2. Picks up the current rack from the support
3. Returns it to its deck position
4. Fetches a fresh rack of the correct volume
5. Loads it onto the support
6. Proceeds with your operation

All of this happens automatically—you just call `transfer_96()` with the appropriate `TrackedTips` instance.

## Setup: Define Both Tip Types

First, set up TrackedTips instances for each volume you'll use:

```python
from pyhamilton import LayoutManager, layout_item, TrackedTips, TipSupportTracker
from pyhamilton.resources import Plate96, Tip96

lmgr = LayoutManager(lay_file)

# Setup 300 µL tips (8 racks on deck)
tracked_tips_300uL = TrackedTips.from_prefix(
    tracker_id="stf_l",
    volume_capacity=300,
    prefix="stf_l",
    count=8,
    tip_type=Tip96,
    lmgr=lmgr
)

# Setup 50 µL tips (8 racks on deck)
tracked_tips_50uL = TrackedTips.from_prefix(
    tracker_id="TIP_50ulF_L",
    volume_capacity=50,
    prefix="TIP_50ulF_L",
    count=8,
    tip_type=Tip96,
    lmgr=lmgr
)
```

Each `TrackedTips` instance manages its own pool of racks. They operate independently but can share the same tip support.

## Setup: Create the Tip Support

The tip support is the physical position where racks are loaded for the 96-channel head to access:

```python
tip_support_resource = layout_item(lmgr, Tip96, 'TipSupport_0001')
tip_support = TipSupportTracker(tip_support_resource)
```

This creates a tracker for the tip support position on your deck.

## Define Your Plates

Set up source and destination plates as usual:

```python
source_plate = layout_item(lmgr, Plate96, 'HSP_Pipette')
dest_plate = layout_item(lmgr, Plate96, 'HSP_Pipette2')
```

## Transfer 1: Using 300 µL Tips

For your first transfer with a larger volume, use the 300 µL tips:

```python
transfer_96(
    ham_int, tracked_tips_300uL, tip_support,
    num_samples=8,
    source_plate=source_plate,
    target_plate=dest_plate,
    volume=100,
    liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty'
)
```

**What happens:**
- The tip support is empty, so a 300 µL rack is automatically loaded
- 1 column (8 tips) is picked up from the rightmost position
- The transfer is performed
- Tips are ejected

## Transfer 2: Switching to 50 µL Tips

For a smaller volume operation, switch to 50 µL tips:

```python
transfer_96(
    ham_int, tracked_tips_50uL, tip_support,
    num_samples=8,
    source_plate=source_plate,
    target_plate=dest_plate,
    volume=20,
    liquid_class='Tip_50ulFilter_Water_DispenseJet_Empty'
)
```

**What happens:**
- The system detects that 300 µL tips are on the support, but 50 µL tips are needed
- The 300 µL rack is automatically picked up and returned to its deck position
- A fresh 50 µL rack is fetched and loaded onto the support
- 1 column of 50 µL tips is picked up
- The transfer is performed
- Tips are ejected

**Important:** Notice you also changed the `liquid_class` parameter to match the 50 µL tips. This is required for proper liquid handling.

## Transfer 3: Switching Back to 300 µL Tips

Switch back to larger tips for another operation:

```python
transfer_96(
    ham_int, tracked_tips_300uL, tip_support,
    num_samples=8,
    source_plate=source_plate,
    target_plate=dest_plate,
    volume=150,
    liquid_class='StandardVolumeFilter_Water_DispenseJet_Empty'
)
```

**What happens:**
- The system detects that 50 µL tips are on the support, but 300 µL tips are needed
- The 50 µL rack is automatically swapped out for a 300 µL rack
- The transfer proceeds with 300 µL tips

## Key Points

**Automatic Detection:**
The `TipSupportTracker` checks the volume capacity of the tips you're requesting against what's currently loaded. When they don't match, it automatically swaps racks.

**No Manual Intervention:**
You never explicitly tell the system to swap racks—it happens transparently based on which `TrackedTips` instance you pass to `transfer_96()`.

**Liquid Class Matching:**
Always ensure your liquid class matches your tip volume:
- 300 µL tips → `'StandardVolumeFilter_Water_DispenseJet_Empty'`
- 50 µL tips → `'Tip_50ulFilter_Water_DispenseJet_Empty'`

**Efficient Partial Rack Usage:**
If you only use a few columns from a rack before switching volumes, those remaining tips aren't wasted—they'll be available next time you switch back to that tip type.

**Multiple Tip Types:**
You can set up as many different tip volumes as needed (50 µL, 300 µL, 1000 µL, etc.) and the system will manage switching between all of them.

