"""Test pyrolyse Pose"""

import pytest
import sys
import importlib

from numpy.testing import assert_array_equal

from pyrosetta import init
from pyrolyse.utils import get_pose


init(set_logging_handler=True)


class TestPDBInfoAttributes:
    from pyrolyse.pdbinfo import (_read_attributes, _read_write_attributes,
                                  _get_set_pairs)

    test_attributes = list(_read_attributes)
    test_write_attr = list(_read_write_attributes)
    test_write_attr += [get_set[0] for get_set in _get_set_pairs]

    pose = get_pose('PYTEST')
    try:
        pdbinfo = pose.pdb_info()
    except TypeError:
        pdbinfo = pose.pdb_info

    # Lazy but gets the work done.
    @pytest.mark.parametrize('attr', test_attributes)
    def test_get_readonly_attributes(self, attr):
        # Can we get every read-only attribute without error?
        getattr(self.pdbinfo, attr)

    @pytest.mark.parametrize('attr', test_attributes)
    def test_set_readonly_attributes(self, attr):
        with pytest.raises(AttributeError):
            setattr(self.pdbinfo, attr, None)

    @pytest.mark.parametrize('attr', test_write_attr)
    def test_get_readwrite_attributes(self, attr):
        getattr(self.pdbinfo, attr)

    @pytest.mark.parametrize('attr', test_write_attr)
    def test_set_readwrite_attr(self, attr):
        old_attr = getattr(self.pdbinfo, attr)
        setattr(self.pdbinfo, attr, old_attr)

    @pytest.mark.parametrize('attr', test_write_attr)
    def test_set_readwrite_attr_wrong_type(self, attr):
        wrong_type = []
        with pytest.raises(TypeError):
            setattr(self.pdbinfo, attr, wrong_type)
