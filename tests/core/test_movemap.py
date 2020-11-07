from pyrolyse import get_pose, init
from pyrosetta.rosetta.protocols.simple_moves import SmallMover


init()


class TestMoveMap:
    from pyrolyse.core.movemap import MoveMap

    mmap = MoveMap(bb=True, chi=True)
    small_mv = SmallMover()
    pose = get_pose('PYTEST')

    def test_init_movemap(self):
        assert self.mmap.get_bb(1)
        assert self.mmap.get_chi(1)
        assert not self.mmap.get_branches(1)
        assert not self.mmap.get_nu(1)
