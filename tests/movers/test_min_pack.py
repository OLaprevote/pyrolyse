from pyrosetta import MoveMap, get_score_function
from pyrolyse import get_pose, init
from pyrolyse.movers.min_pack import MinMover


init()


class TestMinMover():
    pose = get_pose('PYTEST')
    mmap = MoveMap()
    mmap.set_bb(True)
    sfxn = get_score_function()

    min_mv = MinMover(movemap=mmap, sfxn=sfxn)

    def test_attributes(self):
        assert self.min_mv.score_function() == self.min_mv.sfxn
        assert self.min_mv.movemap(self.pose).get_bb(1)
