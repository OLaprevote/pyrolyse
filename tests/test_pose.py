"""Test pyrolyse Pose"""

import pytest
from numpy.testing import assert_array_equal

from pyrosetta import init, get_fa_scorefxn
from pyrolyse.pose import Pose, pose_from_sequence, _read_attributes


init(set_logging_handler=True)

lys_pose = pose_from_sequence('PYTEST')

# Getting membrane_info returns an error as no membrane is set up.
test_attributes = list(_read_attributes)    # _read_attributes is a tuple
test_attributes.remove('membrane_info')


# Check getting and setting monkey-patched attributes.
@pytest.mark.parametrize('attr', test_attributes)
def test_get_readonly_attributes(attr):
    getattr(lys_pose, attr)


@pytest.mark.parametrize('attr', test_attributes)
def test_set_readonly_attributes(attr):
    with pytest.raises(AttributeError):
        setattr(lys_pose, attr, None)


@pytest.mark.parametrize('attr', ['fold_tree', 'pdb_info'])
def test_get_readwrite_attributes(attr):
    getattr(lys_pose, attr)


@pytest.mark.parametrize('attr', ['fold_tree', 'pdb_info'])
def test_set_readwrite_attr(attr):
    old_attr = getattr(lys_pose, attr)
    setattr(lys_pose, attr, old_attr)


@pytest.mark.parametrize('attr', ['fold_tree', 'pdb_info'])
def test_set_readwrite_attr_wrong_type(attr):
    wrong_type = 1
    with pytest.raises(TypeError):
        setattr(lys_pose, attr, wrong_type)


# Check data descriptors.
def test_residues():
    assert lys_pose.residues


def test_reslabels():
    assert lys_pose.reslabels


def test_scores():
    # Killing two birds with one stone by also testing score functions.
    sfxn = get_fa_scorefxn()
    sfxn(lys_pose)
    assert lys_pose.scores


@pytest.mark.parametrize('get_torsion,set_torsion,torsion_list',
                         [('phi', 'set_phi', 'phis'), ('psi', 'set_psi', 'psis')]
                         )
class TestPoseTorsionList:
    lys_pose = pose_from_sequence('PYTEST')

    def test_get(self, get_torsion, set_torsion, torsion_list):
        assert_array_equal(getattr(self.lys_pose, torsion_list), [180.0]*6)

    @pytest.mark.parametrize('resid', list(range(-6, 6)))
    def test_get_item_psis(self, get_torsion, set_torsion, torsion_list, resid):
        """Test __getitem__ of parametrized TorsionList

        Ex: say we test the first residue with attribute phis, the code would
        be equivalent to:
        >>> self.lys_pose.set_psi(1, 0.)
        >>> assert self.lys_pose.psis[0] == self.lys_pose.psi(1)
        """
        if resid >= 0:
            index = resid +1
            value = 0.
        else:
            index = self.lys_pose.size + resid + 1
            value = 33.
        getattr(self.lys_pose, set_torsion)(index, value)
        assert (getattr(self.lys_pose, torsion_list)[resid]
                == getattr(self.lys_pose, get_torsion)(index))

    def test_get_slice_psis(self, get_torsion, set_torsion, torsion_list):
        assert_array_equal(getattr(self.lys_pose, torsion_list)[2:-1],
                           [0.] * 3)

    def test_set_psis(self, get_torsion, set_torsion, torsion_list):
        setattr(self.lys_pose, torsion_list, [60.] * 6)
        for resid in range(1, self.lys_pose.size):
            assert getattr(self.lys_pose, get_torsion)(resid) == 60.

    @pytest.mark.parametrize('resid', list(range(-6, 6)))
    def test_set_item_psis(self, get_torsion, set_torsion, torsion_list, resid):
        if resid >= 0:
            index = resid +1
            value = -30.
        else:
            index = self.lys_pose.size + resid + 1
            value = -60.
        getattr(self.lys_pose, torsion_list)[resid] = value
        assert getattr(self.lys_pose, get_torsion)(index) == value

    def test_set_range_psis(self, get_torsion, set_torsion, torsion_list):
        getattr(self.lys_pose, torsion_list)[2:-1] = [90.] * 3
        for resid in range(3, 5):
            assert getattr(self.lys_pose, get_torsion)(resid) == 90.
        assert_array_equal(getattr(self.lys_pose, torsion_list),
                           [-30., -30., 90., 90., 90., -30.])
