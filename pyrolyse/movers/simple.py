from pyrosetta.rosetta.protocols.simple_moves import SmallMover

from ..pythonizer.mover.simple import (CallableProperty, get_angles_max,
                                       set_angles_max)


__all__ = ['SmallMover']


# TODO find a way to change "path" of class.
# Child class not possible because of virtual functions.

SmallMover.movemap = CallableProperty(SmallMover.movemap.__call__,
                                      SmallMover.movemap.__get__,
                                      SmallMover.movemap)

SmallMover.temperature = property(SmallMover.temperature,
                                  SmallMover.temperature)
SmallMover.nmoves = property(SmallMover.nmoves, SmallMover.nmoves)
SmallMover.angles_max = property(get_angles_max, set_angles_max)

# Make SmallMover callable
SmallMover.__call__ = SmallMover.apply

# TODO wrong doc.
SmallMover.__doc__ = """Mover perturbing phi and psi of a random residue

Such that they create a 'shearing' effect, minimizing the downstream
consequences of this torsional perturbation. The final torsion angle
is subject to a metropolis criterion using the rama score to ensure that
only favorable backbone torsion angles are being selected. The number of
perturbations, and the magnitude of perturbations, and the temperature
in the rama check, can all be modified.

Parameters
----------
1.
movemap_in: pyrosetta.rosetta.core.kinematics.MoveMap
temperature_in: float
nmoves_in: int

2.
arg0: pyrosetta.rosetta.protocols.simple_moves.ShearMover
    If the only input in a ShearMover object, copies it.

Attributes
----------
temperature
nmoves
angles_max: dict

Methods
-------
apply
angle_max
"""
