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
MinMover.__doc__ = """Minimize a pose to a local energy minimum

A protocols::moves::Mover that minimizes a Pose to a local energy minimum by
performing energy minimization of a ScoreFunction over the allowable
degrees of freedom, defined by a MoveMap. The minimization type,
minimization tolerance, and various other options can be also be set.

Modified by pyrolyse at pyrolyse.movers.min_pack

Parameters
----------
1. No parameter
2. str

3.
arg0 : pyrosetta.rosetta.core.kinematics.MoveMap
arg1 : pyrosetta.rosetta.core.scoring.ScoreFunction
arg2 : str
arg3 : float
arg4 : bool

4.
arg0 : pyrosetta.rosetta.core.kinematics.MoveMap
arg1 : pyrosetta.rosetta.core.scoring.ScoreFunction
arg2 : str
arg3 : float
arg4 : bool
arg5 : bool

5.
movemap_in : pyrosetta.rosetta.core.kinematics.MoveMap
scorefxn_in : pyrosetta.rosetta.core.scoring.ScoreFunction
min_type_in : str
tolerance_in : float
use_nb_list_in : bool
deriv_check_in : bool
deriv_check_verbose_in : bool

6. arg0 : pyrosetta.rosetta.protocols.minimization_packing.MinMover

Other parameters
----------------
movemap : pyrosetta.rosetta.core.kinematics.MoveMap, optional
sfxn : pyrosetta.rosetta.core.scoring.ScoreFunction, optional

Attributes
----------
movemap : pyrosetta.rosetta.core.kinematics.MoveMap
sfxn : pyrosetta.rosetta.core.scoring.ScoreFunction

Methods
-------
apply
score_function
    Same as sfxn but as a method.

Examples
--------
>>> min_mv = lys.movers.min_pack.MinMover(movemap=mmap, sfxn=sfxn)
>>> print(min_mv.movemap)
>>> min_mv(pose)
"""
