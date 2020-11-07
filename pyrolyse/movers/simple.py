from pyrosetta.rosetta.protocols.simple_moves import SmallMover, ShearMover

from ..bindings.movers import (_CallableProperty, _get_angles_max,
                                 _set_angles_max)


__all__ = ['SmallMover', 'ShearMover']

# movemap can now also be set using mover.movemap = new_movemap
# One pitt fall, though: writing mover.movemap = pose would
# be the same as writing mover.movemap(pose).
SmallMover.movemap = _CallableProperty(SmallMover.movemap.__call__,
                                      SmallMover.movemap.__get__,
                                      SmallMover.movemap)

SmallMover.temperature = property(SmallMover.temperature, SmallMover.temperature)
SmallMover.nmoves = property(SmallMover.nmoves, SmallMover.nmoves)
SmallMover.angles_max = property(_get_angles_max, _set_angles_max)

SmallMover.__call__ = SmallMover.apply
SmallMover.__doc__ = """Mover randomly perturbing phi and psi of a residue

A mover that makes independent random perturbations of the phi and
psi torsion angles of residue i. It selects residue i at random among
movable residues (set by its MoveMap), and the final torsion angle
is subject to a metropolis criterion using the rama score to ensure that
only favorable backbone torsion angles are being selected. The number of
perturbations, and the magnitude of perturbations, and the temperature
in the rama check, can all be modified.

Parameters
----------
1. No argument
2.
movemap_in: pyrosetta.rosetta.core.kinematics.MoveMap
temperature_in: float
nmoves_in: int

3.
arg0: pyrosetta.rosetta.protocols.simple_moves.SmallMover
    If the only input in a SmallMover object, copies it.

Attributes
----------
temperature
nmoves
angles_max: dict

Methods
-------
apply
angle_max

Examples
--------
>>> mmap = lys.MoveMap(bb=True)
>>> small_mv = lys.movers.simple.SmallMover(mmap, 1., 1)
>>> small_mv.angles_max = {'H':5, 'E': 10, 'L':15}
>>> pose = lys.get_pose('A'*5)
>>> small_mv(pose)
"""


# Make movemap settable
ShearMover.movemap = _CallableProperty(ShearMover.movemap.__call__,
                                      ShearMover.movemap.__get__,
                                      ShearMover.movemap)

ShearMover.temperature = property(ShearMover.temperature, ShearMover.temperature)
ShearMover.nmoves = property(ShearMover.nmoves, ShearMover.nmoves)
ShearMover.angles_max = property(_get_angles_max, _set_angles_max)

# Make ShearMover callable
ShearMover.__call__ = ShearMover.apply

ShearMover.__doc__ = """Mover perturbing torsion angles to create a shear effect

Modified from PyRosetta by pyrolyse.
A mover that perturbs the phi of residue i and the psi of residue
i-1 such that they create a 'shearing' effect, minimizing the downstream
consequences of this torsional perturbation. The final torsion angle
is subject to a metropolis criterion using the rama score to ensure that
only favorable backbone torsion angles are being selected. The number of
perturbations, and the magnitude of perturbations, and the temperature
in the rama check, can all be modified.

Parameters
----------
1. No argument.
2.
movemap_in: pyrosetta.rosetta.core.kinematics.MoveMap
temperature_in: float
nmoves_in: int

3.
arg0: pyrosetta.rosetta.protocols.simple_moves.ShearMover
    If the only input in a ShearMover object, copies it.

Attributes
----------
temperature
nmoves
angles_max: lys.pythonize.movers.simple.AngleMaxDict
    Inherited from dict. Get and set angle_max instead of items.

Methods
-------
apply
angle_max

Examples
--------
>>> mmap = lys.MoveMap(bb=True)
>>> shear_mv = lys.movers.simple.ShearMover(mmap, 1., 1)
>>> shear_mv.angles_max = {'H':5, 'E': 10, 'L':15}
>>> pose = lys.get_pose('A'*5)
>>> shear_mv(pose)
"""
