import pytest
import sys
from pathlib import Path

from pyrolyse.utils import init, get_pose

init()

pdb_path = Path(sys.path[0], 'data', '1ubq.pdb')
@pytest.mark.parametrize('pose_input,kind',
                         [(pdb_path, 'file'), ('PYTEST', 'sequence'),
                          ('1UBQ', 'rcsb')])
def test_get_pose(pose_input, kind):
    assert get_pose(pose_input, kind)

@pytest.mark.parametrize('pose_input', [pdb_path, 'PYTEST', '1UBQ'])
def test_auto_get_pose(pose_input):
    assert get_pose(pose_input)
