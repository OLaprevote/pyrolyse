import pytest

from .simple import SmallMover
from pyrosetta import init, MoveMap, pose_from_sequence, Pose

init()

class TestSmallMover:
    def setup_class(cls):
        mmap = MoveMap()

        cls.pose = pose_from_sequence('LITTLE')

        cls.small = SmallMover(mmap, 1., 1)

    def test_attributes(self):
        assert self.small.temperature == 1.
        assert self.small.nmoves == 1

    def test_set_attributes(self):
        kT = 100.
        n = 2
        new_mm = MoveMap()
        new_mm.set_bb(True)

        self.small.temperature = kT
        assert self.small.temperature == kT

        self.small.nmoves = n
        assert self.small.nmoves == n

        self.small.movemap = new_mm
        output_mmap = self.small.movemap(self.pose)
        assert output_mmap.get_bb(1) is True

    def test_apply_pose(self):
        new_pose = Pose()
        new_pose.assign(self.pose)
        self.small(new_pose)
        assert new_pose.residues[-1].xyz(4) != self.pose.residues[-1].xyz(4)

    @pytest.mark.parametrize('struct', ['H','E','L'])
    def test_get_angles_max(self, struct):
        assert type(self.small.angles_max[struct]) == type(float())

    def test_get_angles_max_indice_error(self):
        with pytest.raises(KeyError):
            self.small.angles_max['X']

    @pytest.mark.parametrize('struct', ['H','E','L'])
    def test_set_angles_max_indice(self, struct):
        angle = 20.
        self.small.angles_max[struct] = angle
        assert self.small.angles_max[struct] == angle

    def test_set_angles_max(self):
        angle = 5.5
        new_angles = {struct: angle for struct in 'HEL'}
        self.small.angles_max = new_angles
        assert self.small.angles_max == new_angles

