from pyhamilton import HamiltonInterface, LayoutManager, layout_item, Tip96, TipType

lmgr = LayoutManager('LSK109_MPH_OxfordNanopore_v0.9.2.lay')
tip_rack = layout_item(lmgr, Tip96, 'TIP_50ulF_L_0001')
tip_support_resource = layout_item(lmgr, Tip96, 'TipSupport_0001')

with HamiltonInterface(windowed=True) as ham_int:
    ham_int.initialize()

    ham_int.set_labware_property(tip_support_resource.layout_name(),'MlStarCore96TipRack', TipType.uL_50)
    # Load tips into the tip rack from the tip support
    ham_int.tip_pick_up_96(tip_support_resource)
    # Drop tips back into the tip support from the tip rack
    ham_int.tip_eject_96(tip_rack)
