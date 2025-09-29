# Precise Reagent Dispensing

NGS protocols usually involve transfering expensive reagents in low (<10uL) volumes
to sample wells. In order to transfer low volumes as precisely and accurately as possible,
there are several techniques we can use:

* Low-volume tips (10 or 50 uL)
* Slow aspiration speed (10-40 uL/s)
* Surface following with cLLD used during aspiration
* Pre-wetting the tip by mixing before aspirating to ensure favorable liquid adhesion 

Here is a pip_transfer function call that uses these:
```python
pip_transfer(
    ham_int, tips=tracked_tips_50uL, source_positions=FastSelect_position,
    dispense_positions=HSP_CPAC_positions, volumes=volumes, prewet_volume=20, prewet_cycles=2,
    liquid_following_aspiration=True, liquid_class='Tip_50ulFilter_Water_DispenseSurface_Empty',
    aspiration_height=0, dispense_height=1
)
```
*See pipette_reagent_to_samples.py for the complete script*