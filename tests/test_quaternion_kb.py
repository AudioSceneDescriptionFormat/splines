# Error case from https://github.com/AudioSceneDescriptionFormat/splines/issues/44

from splines.quaternion import UnitQuaternion, KochanekBartels

def test_create_kb():
    kb = KochanekBartels([
        UnitQuaternion.from_unit_xyzw([0,0,1,0]),
        UnitQuaternion.from_unit_xyzw([0,0,1,0]),
    ])
    assert len(kb.grid) == 2
