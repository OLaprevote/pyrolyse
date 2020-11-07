import pytest

from importlib import reload
from sys import modules

# Avoid using other monkey-patched classes imported in other tests
try:
    reload(modules['pyrosetta'])
except: pass

from pyrosetta import init, MoveMap, Pose, pose_from_sequence
from pyrosetta.rosetta.core.pose import make_pose_from_sequence
from pyrolyse.movers.simple import SmallMover, ShearMover

init(set_logging_handler=True)

class TestSmallMover:
    def setup_class(cls):
        mmap = MoveMap()

        cls.pose = Pose()
        make_pose_from_sequence(cls.pose, 'PYTEST', "fa_standard", True)

        cls.small_mv = SmallMover(mmap, 1., 1)

    def test_attributes(self):
        assert self.small_mv.temperature == 1.
        assert self.small_mv.nmoves == 1

    def test_set_attributes(self):
        kT = 100.
        n = 2
        new_mm = MoveMap()
        new_mm.set_bb(True)

        self.small_mv.temperature = kT
        assert self.small_mv.temperature == kT

        self.small_mv.nmoves = n
        assert self.small_mv.nmoves == n

        self.small_mv.movemap = new_mm
        output_mmap = self.small_mv.movemap(self.pose)
        assert output_mmap.get_bb(1) is True

    def test_apply_pose(self):
        new_pose = Pose()
        new_pose.assign(self.pose)
        self.small_mv(new_pose)
        assert new_pose.residue(6).xyz(4) != self.pose.residue(6).xyz(4)

    @pytest.mark.parametrize('struct', ['H','E','L'])
    def test_get_angles_max(self, struct):
        assert type(self.small_mv.angles_max[struct]) == type(float())

    def test_get_angles_max_indice_error(self):
        with pytest.raises(KeyError):
            self.small_mv.angles_max['X']

    @pytest.mark.parametrize('struct', ['H','E','L'])
    def test_set_angles_max_indice(self, struct):
        angle = 20.
        self.small_mv.angles_max[struct] = angle
        assert self.small_mv.angles_max[struct] == angle

    def test_set_angles_max(self):
        angle = 5.5
        new_angles = {struct: angle for struct in 'HEL'}
        self.small_mv.angles_max = new_angles
        assert self.small_mv.angles_max == new_angles


class TestShearMover:
    def setup_class(cls):
        mmap = MoveMap()

        cls.pose = Pose()
        make_pose_from_sequence(cls.pose, 'PYTEST', "fa_standard", True)

        cls.shear_mv = ShearMover(mmap, 1., 1)

    def test_attributes(self):
        assert self.shear_mv.temperature == 1.
        assert self.shear_mv.nmoves == 1

    def test_set_attributes(self):
        kT = 100.
        n = 2
        new_mm = MoveMap()
        new_mm.set_bb(True)

        self.shear_mv.temperature = kT
        assert self.shear_mv.temperature == kT

        self.shear_mv.nmoves = n
        assert self.shear_mv.nmoves == n

        self.shear_mv.movemap = new_mm
        output_mmap = self.shear_mv.movemap(self.pose)
        assert output_mmap.get_bb(1) is True

    def test_apply_pose(self):
        new_pose = Pose()
        new_pose.assign(self.pose)
        self.shear_mv(new_pose)
        assert new_pose.residue(6).xyz(4) != self.pose.residue(6).xyz(4)

    @pytest.mark.parametrize('struct', ['H','E','L'])
    def test_get_angles_max(self, struct):
        assert type(self.shear_mv.angles_max[struct]) == type(float())

    def test_get_angles_max_indice_error(self):
        with pytest.raises(KeyError):
            self.shear_mv.angles_max['X']

    @pytest.mark.parametrize('struct', ['H','E','L'])
    def test_set_angles_max_indice(self, struct):
        angle = 20.
        self.shear_mv.angles_max[struct] = angle
        assert self.shear_mv.angles_max[struct] == angle

    def test_set_angles_max(self):
        angle = 5.5
        new_angles = {struct: angle for struct in 'HEL'}
        self.shear_mv.angles_max = new_angles
        assert self.shear_mv.angles_max == new_angles

