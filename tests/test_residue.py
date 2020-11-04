"""Test pyrolyse Pose"""

import pytest
import sys
import importlib

from numpy.testing import assert_array_equal

from pyrosetta import init, pose_from_sequence


init(set_logging_handler=True)


class TestResidueAttributes:
    from pyrolyse.residue import _read_attributes, _read_write_attributes
    from pyrolyse.pose import pose_from_sequence

    test_attributes = list(_read_attributes)
    test_attributes.remove('carbohydrate_info')
    test_write_attr = list(_read_write_attributes)
    ros_pose = pose_from_sequence('PYTEST')
    lys_res = ros_pose.residue(2)

    # Lazy but gets the work done.
    @pytest.mark.parametrize('attr', test_attributes)
    def test_get_readonly_attributes(self, attr):
        # Can we get every read-only attribute without error?
        print(attr)
        getattr(self.lys_res, attr)

    @pytest.mark.parametrize('attr', test_attributes)
    def test_set_readonly_attributes(self, attr):
        with pytest.raises(AttributeError):
            setattr(self.lys_res, attr, None)

    @pytest.mark.parametrize('attr', test_write_attr)
    def test_get_readwrite_attributes(self, attr):
        getattr(self.lys_res, attr)

    @pytest.mark.parametrize('attr', test_write_attr)
    def test_set_readwrite_attr(self, attr):
        old_attr = getattr(self.lys_res, attr)
        setattr(self.lys_res, attr, old_attr)

    @pytest.mark.parametrize('attr', test_write_attr)
    def test_set_readwrite_attr_wrong_type(self, attr):
        wrong_type = None
        with pytest.raises(TypeError):
            setattr(self.lys_res, attr, wrong_type)
