"""
TODO

Torsion list
------------
alpha
beta
delta
epsilon
gamma
mu
omega
theta
zeta
phi
psi

Read/Set property
-----------------
fold_tree   # Read & set property
pdb_info    # Read & set property

Other ideas
-----------

conformation # Read-only attribute
set_new_conformation # Setter ???

append_residue # Summarize all append_residue methods?
batch_get_xyz and batch_set_xyz # See how it works
insert_residue  # See if possible to summarize insert_residue_by_bond, insert_residue_by_jump
                # and insert_residue_by_atoms

chi     # (chino: int, seqpos: int) a numpy array? -> Wouldn't work, not same dimension.

aa # Transform into a list with callables.


pdb_rsd     # Numpy-ish (takes a tuple (chain, resNo) as argument)

replace_residue     # To use as a residue setter?

secstruct   # pose.secstruct[24] returns secstruct of residue 24?
set_secstruct   # then settable for range of residues.
"""

from pyrosetta.rosetta.core.pose import (make_pose_from_sequence, PDBInfo,
                                         Pose)

import numpy as np

__all__ = ['Pose', 'pose_from_sequence']


class TorsionList(list):
    def __init__(self, instance, torsion_setter, *args):
        self.instance = instance
        self.torsion_setter = torsion_setter
        super().__init__(*args)

    def __setitem__(self, index, value):
        if isinstance(index, int):
            self.torsion_setter(self.instance, index+1, value)

        # Slices are obtained when using semi-colons to get items, e.g. a[2:8:3]
        elif isinstance(index, slice):
            # Find concerned residue indexes
            residues = range(1, self.instance.size+1)[index]
            for resid, new_psi in zip(residues, value):
                self.torsion_setter(self.instance, resid, new_psi)

        else: raise TypeError(("TorsionList indices must be integers or"
                               " slices, not {}").format(type(index).__name__))


class PsiList(list):
    def __init__(self, instance, *args):
        self.instance = instance
        super().__init__(*args)

    def __setitem__(self, index, value):
        if isinstance(index, int):
            self.instance.set_psi(index+1, value)

        # Slices are obtained when using semi-colons to get items, e.g. a[2:8:3]
        elif isinstance(index, slice):
            # Find concerned residue indexes
            residues = range(1, self.instance.size+1)[index]
            for resid, new_psi in zip(residues, value):
                self.instance.set_psi(resid, new_psi)

        else: raise TypeError(("TorsionList indices must be integers or"
                               " slices, not {}").format(type(index).__name__))


def get_psis(instance):
    psi_list = [instance.psi(resid) for resid in range(1, instance.size+1)]
    return TorsionList(instance, Pose.set_psi, psi_list)


def set_psis(instance, new_psis):
    for resid in range(instance.size):
        instance.set_psi(resid+1, new_psis[resid])


def pose_from_sequence(seq, res_type="fa_standard", auto_termini=True):
    """Returns a pose from single-letter protein sequence. /!\ From pyrolyse

    Returns a Pose object generated from a single-letter sequence of amino acid
    residues in <seq> using the <res_type> ResidueType and creates N- and C-
    termini if <auto_termini> is set to True.

    Unlike make_pose_from_sequence(), this method generates a default PDBInfo
    and sets all torsion angles to 180 degrees.

    Example
    -------
    pose = pose_from_sequence("THANKSEVAN")

    See also
    --------
    Pose
    make_pose_from_sequence()
    pose_from_file()
    pose_from_rcsb()
    """
    # Needed corrections to work with pyrolyse Pose
    pose = Pose()
    make_pose_from_sequence(pose, seq, res_type, auto_termini)

    for i in range(0, pose.total_residue):
        res = pose.residue(i + 1)
        if not res.is_protein() or res.is_peptoid() or res.is_carbohydrate():
            continue

        pose.set_phi(i + 1, 180)
        pose.set_psi(i + 1, 180)
        pose.set_omega(i + 1, 180)

    # Empty PDBInfo (rosetta.core.pose.PDBInfo()) is not correct here;
    # we have to reserve space for atoms....
    pose.pdb_info(PDBInfo(pose))
    pose.pdb_info().name(seq[:8])

    return pose

# Read-only attributes
Pose.const_data_cache = property(Pose.const_data_cache)
Pose.data = property(Pose.data)
Pose.dof = property(Pose.dof)
Pose.energies = property(Pose.energies)
Pose.atom_tree = property(Pose.atom_tree)
Pose.membrane_info = property(Pose.membrane_info)
Pose.num_chains = property(Pose.num_chains)
Pose.num_jump = property(Pose.num_jump)
Pose.observer_cache = property(Pose.observer_cache)
Pose.reference_pose_set = property(Pose.reference_pose_set)
Pose.reference_pose_set_cop = property(Pose.reference_pose_set_cop)
Pose.sequence = property(Pose.sequence)
Pose.size = property(Pose.size)
Pose.total_atoms = property(Pose.total_atoms)
Pose.total_residue = property(Pose.total_residue)

# Torsion lists
Pose.psis = property(get_psis, set_psis)
