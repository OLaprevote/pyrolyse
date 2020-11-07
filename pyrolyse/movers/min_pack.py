from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from ..bindings.movers import _CallableProperty


__all__ = ['MinMover']

# movemap can now also be set using mover.movemap = new_movemap
MinMover.movemap = _CallableProperty(MinMover.movemap.__call__,
                                      MinMover.movemap.__get__,
                                      MinMover.movemap)
MinMover.sfxn = property(MinMover.score_function, MinMover.score_function)

_old_init_minmover = MinMover.__init__


# __init__ of MinMover permits to start it with score_function and
# movemap, but requires a load of other arguments I don't think I care
# about. However this way one can optionally already enter movemap and
# sfxn.
# Plus: it's not called movemap_in and scorefxn_in anymore.
def _init_minmover(self, *args, movemap=None, sfxn=None,**kwargs):
    _old_init_minmover(self, *args, **kwargs)
    self.movemap = movemap
    self.sfxn = sfxn


MinMover.__init__ = _init_minmover
