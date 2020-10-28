"""Test pyrolyse Pose"""

import pytest
from numpy.testing import assert_array_equal

from pyrosetta import init
from pyrolyse.pose import Pose, pose_from_sequence

init(set_logging_handler=True)


class TestPose:
    lys_pose = pose_from_sequence('PYTEST')

    def test_get_psis(self):
        assert_array_equal(self.lys_pose.psis, [180.0] * 6)

    def test_get_item_psis(self):
        self.lys_pose.set_psi(1, 0)
        assert self.lys_pose.psis[0] == self.lys_pose.psi(1)

    def test_get_range_psis(self):
        assert_array_equal(self.lys_pose.psis[2:5], [180.0] * 3)

    def test_set_psis(self):
        self.lys_pose.psis = [60.] * 6
        for resid in range(1, self.lys_pose.size):
            assert self.lys_pose.psi(resid) == 60.

    def test_set_item_psis(self):
        self.lys_pose.psis[0] = -30.
        assert self.lys_pose.psi(1) == -30.
        assert_array_equal(self.lys_pose.psis, [-30., 60., 60., 60., 60., 60.])

    def test_set_range_psis(self):
        self.lys_pose.psis[1:4] = [90.] * 3
        for resid in range(2, 5):
            assert self.lys_pose.psi(resid) == 90.
        assert_array_equal(self.lys_pose.psis, [-30., 90., 90., 90., 60., 60.])
