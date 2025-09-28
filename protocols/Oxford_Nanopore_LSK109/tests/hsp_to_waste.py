from pyhamilton import HamiltonInterface, LayoutManager, Plate96, layout_item
from pyhamilton.transport import transport_resource, GripDirection, GrippedResource

lmgr = LayoutManager('LSK109_MPH_OxfordNanopore_v0.9.2.lay')

HSP_Pipette2 = layout_item(lmgr, Plate96, 'HSP_Pipette2')  
HSP_Waste = layout_item(lmgr, Plate96, 'HSP_Waste')

with HamiltonInterface(windowed=True) as ham_int:
    ham_int.initialize()
    transport_resource(ham_int, HSP_Pipette2, HSP_Waste,
                       core_gripper=True, resource_type=GrippedResource.PCR)