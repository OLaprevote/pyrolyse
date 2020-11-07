import pytest

from importlib import reload
from sys import modules

# Try to avoid using other monkey-patched classes.
try:
    reload(modules['pyrosetta'])
except: pass

from pyrosetta import init, MoveMap, Pose
from pyrosetta.rosetta.core.pose import make_pose_from_sequence
from pyrosetta.rosetta.protocols.simple_moves import SmallMover

from pyrolyse.movers import SequenceMover, RepeatMover


init(set_logging_handler=True)


class TestSequenceMover:
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


class TestRepeatMover:
    def setup_class(cls):
        mmap = MoveMap()
        mmap.set_bb(True)

        cls.pose = Pose()
        make_pose_from_sequence(cls.pose, 'PYTEST', "fa_standard", True)

        cls.small_mv = SmallMover(mmap, 1., 1)
        cls.rep_mv = RepeatMover(cls.small_mv, 3)

    def test_apply_pose(self):
        new_pose = self.pose.clone()
        self.rep_mv(new_pose)
        assert new_pose.residue(6).xyz(4) != self.pose.residue(6).xyz(4)
