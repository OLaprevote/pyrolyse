from pyrosetta.rosetta.protocols.moves import SequenceMover, RepeatMover


__all__ = ('SequenceMover', 'RepeatMover')

_old_init_seq_mover = SequenceMover.__init__

def _init_seq_mover(self, *args, movers=(), **kwargs):
    _old_init_seq_mover(self, *args, **kwargs)
    for mover in movers:
        try:
            self.add_mover(mover)
        except TypeError:
            # In case movers are weighted.
            try:
                self.add_mover(*mover)
            except TypeError:
                self.add_mover(**mover)


SequenceMover.__init__ = _init_seq_mover
SequenceMover.__call__ = SequenceMover.apply
SequenceMover.__doc__ = """A Mover that iterates through a vector of Movers

Modified from PyRosetta by pyrolyse. Applies each mover sequentially.
Instances are callable, equivalent to using apply method.

Parameters
----------
1.
No argument

2.
ms: bool

3.
mover1: pyrosetta.rosetta.protocols.moves.Mover
mover2: pyrosetta.rosetta.protocols.moves.Mover

4.
mover1: pyrosetta.rosetta.protocols.moves.Mover
mover2: pyrosetta.rosetta.protocols.moves.Mover
mover3: pyrosetta.rosetta.protocols.moves.Mover

5.
arg0: pyrosetta.rosetta.protocols,moves.SequenceMover

movers: iterable object, optional
    Sequence of movers, arguments or kwargs added sequentially to the
    instance.
    Ex: (Mover(), (Mover(), 2.), {'mover_in': Mover(), 'weight_in': 2.})
Examples
--------
>>> mmap = lys.MoveMap(bb=True)
>>> small_mov = lys.movers.simple.SmallMover(mmap, 1., 1)   # Random phi/psi change
>>> pose = lys.get_pose('A'*5)
>>> two_small_mov = lys.SequenceMover([small_mov, small_mov])
>>> two_small_mov(pose)    # Randomly modifies 2 times phis and psis of pose
"""


RepeatMover.__call__ = RepeatMover.apply
RepeatMover.__doc__ = """Mover to repeat an input Mover a specified number of times

Modified from PyRosetta by pyrolyse. Instances are callable, equivalent
to using apply method.

Parameters
----------
1. Nothing
2.
mover_in: pyrosetta.rosetta.protocols.moves.Mover
    Mover to repeat
nmoves_in: int
    Number of times to repeat
3.
arg0: pyrosetta.rosetta.protocols.moves.RepeatMover
    RepeatMover to copy

Examples
--------
>>> mmap = lys.MoveMap(bb=True)
>>> small_mov = lys.movers.simple.SmallMover(mmap, 1., 1) # Random phi/psi change
>>> pose = lys.get_pose('A'*5)
>>> five_small_mov = lys.RepeatMover(small_mov, 5)
>>> five_small_mov(pose)  # Randomly modifies 5 times phis and psis of pose
"""
