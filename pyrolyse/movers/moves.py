from pyrosetta.rosetta.protocols.moves import SequenceMover


__all__ = ('SequenceMover')

_old_init_seq_mover = SequenceMover.__init__

def _init_seq_mover(self, *args, movers=()):
    _old_init_seq_mover(self, *args)
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
"""
