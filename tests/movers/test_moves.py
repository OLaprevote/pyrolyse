import pytest

from importlib import reload
from sys import modules

# Try to avoid using other monkey-patched classes imported in other tests
# "It is not very efficient"
try:
    reload(modules['pyrosetta'])
except: pass

from pyrosetta import init, MoveMap, Pose
from pyrosetta.rosetta.core.pose import make_pose_from_sequence
from pyrosetta.rosetta.protocols.simple_moves import SmallMover

from pyrolyse.movers import SequenceMover

init(set_logging_handler=True)

class TestSmallMover:

    def setup_class(cls):
        mmap = MoveMap()
        mmap.set_bb(True)

        cls.pose = Pose()
        make_pose_from_sequence(cls.pose, 'PYTEST', "fa_standard", True)

        cls.small = SmallMover(mmap, 1., 1)
        cls.seqmov = SequenceMover()

    def test_sequence_movers(self):
        self.seqmov.assign(SequenceMover(movers=[self.small]*3))
        assert self.seqmov.movers()[0]==self.small

    def test_apply_pose(self):
        new_pose = self.pose.clone()
        self.seqmov(new_pose)
        assert new_pose.residue(6).xyz(4) != self.pose.residue(6).xyz(4)
