from pyrosetta.rosetta.core.kinematics import MoveMap


_old_init_movemap = MoveMap.__init__


def _init_movemap(self, *args, bb=False, branches=False, chi=False, jump=False,
              nu=False, **kwargs):
    _old_init_movemap(self, *args, **kwargs)
    self.set_bb(bb)
    self.set_branches(branches)
    self.set_chi(chi)
    self.set_jump(jump)
    self.set_nu(nu)


MoveMap.__init__ = _init_movemap
MoveMap.__doc__ = """Class specifying DOFs to be flexible or fixed

Currently there are two groups of data, one is a residue-based Torsion
definition, such as BB, CHI, NU, and JUMP; the other is an atom-based DOF
definition, such as bond length D, bond angle THETA, and torsion angle PHI,
which are used in the AtomTree. MoveMap does not automatically handle
conversion from one group to the other, i.e., setting PHI false for
DOF_type does not affect setting for BB and CHI torsion though they are
PHIs in atom-tree.

Within each group, there are multiple levels of control
(from general/high to specific/lower):


Torsion-based: TorsionType(BB, CHI, NU, BRANCH, JUMP) -> MoveMapTorsionID
(BB, CHI of one residue) -> TorsionID ( BB torsion 2 or CHI torsion 3 of
one residue)


DOF-base: DOF_type( D, THETA, PHI ) -> DOF_ID (D, THETA, PHI of one atom)

Settings for each level are stored in a map structure and they are only
added to the map when setting methods are invoked. As a result, MoveMap
does not behave like a "Boolean vector", which always contains setting for
each residue or atom in a conformation. Setting for a higher level will
override setting for lower levels (remove it from map); Similarly, when
querying a lower level finds no setting, it will check setting for its
higher level. For example, setting TorsionType BB to be true will remove
any data of BB setting for a residue or a specific BB torsion (such as
backbone psi) in a residue. And querying the setting for BB torsion 2 of
residue 4 will first check if there is any specific setting, if not, it will
check if there is a setting for all BB torsions for residue 4, if not
again, it will use the setting for BB torsions for all residues.

Modified by pyrolyse at pyrolise.core.movemap.

Parameters
----------
1. No argument
2. arg0 : pyrosetta.rosetta.core.kinematics.MoveMap

Other Parameters
----------------
bb : bool, optional
    Set whether or not BB TorsionTypes are movable
branches : bool, optional
    Set whether or not BRANCH TorsionTypes are movable
chi : bool, optional
    Set whether or not CHI TorsionTypes are movable
jump : bool, optional
    Set whether or not BRANCH TorsionType is movable
nu : bool, optional
    Set whether or not NU TorsionTypes are movable

See also
--------
Pose
MinMover
ShearMover
SmallMover

Examples
--------
>>> mmap = lys.MoveMap()
>>> mmap.set_bb(True)
>>> mmap2 = lys.MoveMap(bb=True)
>>> mmap3 = lys.MoveMap(mmap2, chi=True, jump=True)
"""
