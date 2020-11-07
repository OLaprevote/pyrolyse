from pyrosetta.rosetta.protocols.minimization_packing import MinMover


MinMover.movemap = property(MinMover.movemap, MinMover.movemap)
MinMover.sfxn = property(MinMover.score_function, MinMover.score_function)

_old_init_minmover = MinMover.__init__


def _init_minmover(self, *args, movemap=None, sfxn=None,**kwargs):
    _old_init_minmover(self, *args, **kwargs)
    self.movemap = movemap
    self.sfxn = sfxn


MinMover.__init__ = _init_minmover
