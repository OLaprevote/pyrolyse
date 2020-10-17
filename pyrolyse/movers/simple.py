import copy
from pyrosetta.rosetta.protocols.simple_moves import SmallMover


class _AngleDict(type(dict())):
    """Dict with setitem setting max_angle values from an object

    Parameters:
    -----------
    obj
        Object with an angle_max method.
    """
    def __init__(self, obj):
        self.obj = obj
        angles_dict = {sec: obj.get_angle_max(sec) for sec in 'HEL'}
        type(dict()).__init__(self, **angles_dict)

    def __setitem__(self, index, value):
        self.obj.angle_max(index, value)


def _get_angles_max(self):
    return(_AngleDict(self))


def _set_angles_max(self, angles):
    for sec in angles.keys():
        self.angle_max(sec, angles[sec])

# TODO
# Have movemap as a settable function, so that
# >>> smallmover.movemap(pose) returns MoveMap object
# and the movemap can be set like
# >>> smallmover.movemap = MoveMap()
#
# def _get_movemap(self):
#     return SmallMover.movemap
#
#
# def _set_movemap(self, value):
#     super().movemap(value)


# TODO find a way to change "path" of class.
# Child class not possible because of virtual functions.
SmallMover.movemap = property(_get_movemap, _set_movemap)
SmallMover.temperature = property(SmallMover.temperature, SmallMover.temperature)
SmallMover.nmoves = property(SmallMover.nmoves, SmallMover.nmoves)
SmallMover.angles_max = property(_get_angles_max, _set_angles_max)

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
