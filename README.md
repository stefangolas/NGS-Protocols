# Python for NGS Automation

This library provides 5 complete automated sequencing prep protocols for the Hamilton NGS platform written in PyHamilton.
These protocols make extensive use of new PyHamilton features specifically designed to allow for easy development of NGS protocols.
Here are some of the main features that have been added to support these protocols:

-   Tip tracking
-   Reagent consumption tracking
-   User prompts for protocol step selection and deck loading
-   JSON-defined liquid class import tool
-   In-method liquid class parameter updates for iterative liquid class development
-   Type checking for command parameters
-   Flexible liquid handling wrapper functions for precise transfers
-   TADM curve report generation
-   Partial rack tip-pickups with the MPH (96-channel head) using column offsets
-   Transport controller with pre-defined gripper parameters for different resources
-   Expanded device library including ODTC, CPAC, and heater shakers


The five protocols provided in this library are:
-   10X GEM-X Single Cell 3' sequencing prep
-   PacBio HiFi Plex
-   Oxford Nanopore LSK109/LSK114
-   QIAseq RNA Fusion
-   KAPA Roche HyperPrep HyperPlus

These protocols are intended to reproduce the exact functionality of validated methods written in
Venus, but they were not directly tested using biological samples. They are intended to be illustrative
guides and starting points for writing NGS sequencing prep protocols in PyHamilton, rather than
fully functional protocols that work out-of-the-box. Please make sure to fully understand
and test any code in this library before deploying it in production.

## Installation

To use the scripts in this library

1. Install the latest version of PyHamilton.
2. Run the asset import script with `install_assets.py`. Note that this will copy files from `/assets` into your `Hamilton` directory. Feel free to do this manually instead.


## PyHamilton NGS Features

This library contains examples of many new PyHamilton features. Here are explanations and simple tutorials for how to use them.

### Tip tracking
The new `TrackedTips` class lets us define a set of tip racks that we can draw from incrementally with `tracked_tip_pick_up`, without having
to manually specify tip positions. We also persistently record the state of our tracked tips in a SQLite database.

```python
tracked_tips_50uL = TrackedTips.from_prefix(
    tracker_id="TIP_50uLF_L",
    volume_capacity=50,
    prefix="TIP_50uLF_L",
    count=8,
    tip_type=Tip96,
    lmgr=self.lmgr)

tracked_tip_pick_up(ham_int, tip_tracker, num_tips = 8)
```

### Reagent consumption tracking
We can use members of the `TrackedReagentVessel` class to specify the contents of a container, record the total amount of
volume aspirated during a run, and generate JSON reports about reagent consumption. Note that `tracked_reagent_aspirate`
must be used in order for the volume of an aspiration command to be tracked.

```python
# First we define a container using a ReagentTracked... class
MagBeads_Container = layout_item(lmgr, ReagentTrackedReservoir60mL, 'rgt_cont_60ml_BC_A00_0001')

# Next we assign a reagent to indexes within that container via assign_reagent_map, which returns a positions list ready to be passed to a command
magbead_positions = MagBeads_Container.assign_reagent_map('MagBeads', range(8))

# Here's an example with a plate instead of a trough
CPAC_Reagent_Plate = layout_item(lmgr, ReagentTrackedPlate96, 'CPAC_HSP_0001')
ER_Mix_positions = CPAC_Reagent_Plate.assign_reagent_map('ER_Mix', [0])
```

### User interface for running protocols
Protocols can be run via a step-selection interface that allows for simulated testing, reagent consumption calculations, and deck loading dialogues.

[See this guide for more details.](/guides/protocol_gui)

### Flexible liquid handling wrapper functions for precise transfers
We now have several wrappers for frequently used transfer operations. Refer to the scripts in the `examples/pipetting` folder
for a thorough exploration of use-cases and parameters.

*Transfering from a reagent source to a partial plate of samples*
```python
lmgr = LayoutManager('PacBio_MultiPlexLibraryPrepDeck_v1.2.lay')
MagBeads_Container = layout_item(lmgr, ReagentTrackedReservoir60mL, 'rgt_cont_60ml_BC_A00_0001')
MIDI_OffMagnet = layout_item(lmgr, Plate96, 'MIDI_Pipette')

tracked_tips_300uL = TrackedTips.from_prefix(
            tracker_id="STF_L",
            volume_capacity=300,
            prefix="STF_L",
            count=8,
            tip_type=Tip96,
            lmgr=lmgr)

magbead_positions =  MagBeads_Container.assign_reagent_map('MagBeads', range(8))
sample_poss = [(MIDI_OffMagnet, i) for i in range(num_samples)]
vols = [magbead_volume] * num_samples

# The source position is repeated for every aspirate, while the dispense position is incremented
pip_transfer(ham_int, tracked_tips_300uL, magbead_positions, sample_poss, vols, dispense_height=1,
                         prewet_cycles = 1,
                         prewet_volume = 20,
                         liquid_class='PacBio_SVF_ProNexPurificationBeads_SurfaceEmpty_v3')
```

*Aspirating supernatant from beads*
```python
# Double aspiration is used to allow residual supernatant to settle at the bottom of the container.
double_aspirate_supernatant_96(ham_int, tracked_tips_300uL, tip_support, num_samples, MIDI_OnMagnet, 
                                          LiquidWaste, first_volume=270, second_volume=30, 
                                          second_aspiration_height=0, liquid_class='StandardVolumeFilter_Water_DispenseSurface_Empty')
```

### Partial rack tip-pickups with the MPH (96-channel head) using column offsets
We often use the 96-channel head to transfer a sample mixture from one plate to another, such as during a bead clean-up.
In order to economize on the use of tips for transfering a partially filled plate, we can load tips onto the Tip Support
holder and pick up a partial rack of tips with an offset based on the number of columns we are picking up.

There is a lot of code under the hood for supporting this, but from a user perspective you just need to instantiate
the appropriate TipTracker, the TipSupportTracker, and call the `mph_tip_pickup_support` command with the number
of samples you are transfering.

```python
tips = TrackedTips.from_prefix(
            tracker_id="STF_L",
            volume_capacity=300,
            prefix="STF_L",
            count=8,
            tip_type=Tip96,
            lmgr=lmgr)

tip_support_resource = layout_item(lmgr, Tip96, 'TipSupport_0001')
tip_support = TipSupportTracker(tip_support_resource)

mph_tip_pickup_support(ham_int, tips, tip_support, num_tips=24)
```

The TipSupportTracker automates loading in a new rack to meet the desired quantity of tips,
and switching between racks of different types when necessary.

### Transport controller
The `transport_resource` class provides a concise and generic wrapper function for moving plates and lids, with specific parameters
for different use-contexts applied under-the-hood. The user specifies a GrippedResource type, whether they are using the iswap or core gripper,
and if using the iswap, also specifies a grip direction. This wrapper manages the complexity and details of transport operations, allowing for highly readable
transports to be specified.

*core gripper transport to magnet stand*
```python
transport_resource(ham_int, MIDI_OffMagnet, MIDI_OnMagnet, 
                             resource_type=GrippedResource.MIDI, core_gripper=True)
```

*iswap transport to thermal cycler*
```python
transport_resource(ham_int, HSP_Pipette, HSP_ODTC, 
                    grip_direction=GripDirection.RIGHT, resource_type=GrippedResource.PCR, iswap=True)
```

### JSON liquid class imports
Instead of having to define liquid classes strictly through the Venus Liquid Editor, you can
specify these in JSON and import them into your liquid class database through a PyHamilton command.

*my_liquid_class.json*
```json
    [
        {
        "name": "Beads",
        "aspirate": {
            "FLOW_RATE": 100,
            "MIX_FLOW_RATE": 100,
            "AIR_TRANSPORT_VOLUME": 5,
            "BLOW_OUT_VOLUME": 30,
            "SWAP_SPEED": 2,
            "SETTLING_TIME": 1,
            "OVER_ASPIRATE_VOLUME": 0,
            "CLOT_RETRACT_HEIGHT": 0
        },
        "dispense": {
            "FLOW_RATE": 180,
            "MIX_FLOW_RATE": 1,
            "AIR_TRANSPORT_VOLUME": 5,
            "BLOW_OUT_VOLUME": 30,
            "SWAP_SPEED": 1,
            "SETTLING_TIME": 0,
            "STOP_FLOW_RATE": 100,
            "STOP_BACK_VOLUME": 0
        },
        "tip_type": {
            "volume": 300,
            "has_filter": true
        },
        "dispense_mode": "Jet Empty",
        "correction_curve": {
            "nominal": [10, 50, 100, 200, 300],
            "corrected": [9.8, 49.5, 99.2, 198.5, 298.8]
        }
    }
]
```

*import_lc.py*
```python
from pyhamilton import HamiltonInterface, copy_liquid_class, create_liquid_class_from_json

with HamiltonInterface(windowed=True) as ham_int:
    ham_int.initialize()
    create_liquid_class_from_json(ham_int, 'my_lc.json')
```

### In-method liquid class parameter updates
We can make in-method updates to liquid class parameters. This is exceptionally useful for custom scripts for
iterative liquid class development. Under the hood, this is also how the above import function works.

```python
from pyhamilton import HamiltonInterface, copy_liquid_class, update_aspirate_parameter, AspirateParameter

# Pick a liquid class to modify, and copy that instead of modifying directly.
copy_liquid_class(ham_int, 'template_liquid_class', 'new_liquid_class')

# Make sure to call the relevant AspirateParameter or DispenseParameter attribute
update_aspirate_parameter(ham_int, 'new_liquid_class', AspirateParameter.FLOW_RATE, 20)
```
