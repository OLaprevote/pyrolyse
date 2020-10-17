import copy
from pyrosetta.rosetta.protocols.simple_moves import SmallMover as _SmallMover

SmallMover = copy.deepcopy(_SmallMover)

class _AngleDict(type(dict())):
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


def _get_movemap(self):
    return SmallMover.movemap


def _set_movemap(self, value):
    super().movemap(value)


# TODO find a way to change "path" of class.
# Creating a child class of SmallMover crashes because of C++ virtual functions.
# See for yourself:
# class Small(SmallMover):
#    pass

# Add attributes to SmallMover.
# SmallMover.movemap = property(SmallMover.movemap, SmallMover.movemap)
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
movemap
temperature
nmoves
angles_max: dict

Methods
-------
apply
angle_max
"""

if __name__ == '__main__':
    import pyrosetta as ros

    mm = ros.MoveMap()
    mover = SmallMover(mm, 1., 1)
    mover2 = _SmallMover(mm, 1., 1)
    print(mover.temperature, mover2.temperature)
