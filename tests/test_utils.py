import pytest

from pyrolyse.utils import init, get_pose

init()


@pytest.mark.parametrize('pose_input,kind',
                         [('data/1ubq.pdb', 'file'), ('PYTEST', 'sequence'),
                          ('1UBQ', 'rcsb')])
def test_get_pose(pose_input, kind):
    assert get_pose(pose_input, kind)

@pytest.mark.parametrize('pose_input', ['data/1ubq.pdb', 'PYTEST', '1UBQ'])
def test_auto_get_pose(pose_input):
    assert get_pose(pose_input)
