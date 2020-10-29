"""
TODO

Other ideas
-----------
conformation # Read-only attribute
set_new_conformation # Setter?

Check residues, reslabels and scores.

append_residue # Summarize all append_residue methods?
batch_get_xyz and batch_set_xyz # See how it works
insert_residue  # See if possible to summarize insert_residue_by_bond, insert_residue_by_jump
                # and insert_residue_by_atoms

chi     # (chino: int, seqpos: int) a numpy array? -> Wouldn't work,
        # not same dimension.
        # Could create a list of list with tuples and tuples of slices
        # understanding.

aa # Transform into a list with callables.

pdb_rsd     # Numpy-ish (takes a tuple (chain, resNo) as argument)

replace_residue     # To use as a residue setter?

secstruct   # pose.secstruct[24] returns secstruct of residue 24?
set_secstruct   # then settable for range of residues.
"""

from pyrosetta.rosetta.core.pose import (make_pose_from_sequence,
                                         PDBInfo, Pose)

from .pythonizer.pose import TorsionList, torsion_list_property


__all__ = ['Pose', 'pose_from_sequence']


def pose_from_sequence(seq, res_type="fa_standard", auto_termini=True):
    """Returns a pose from single-letter protein sequence.

    Modified from pyrosetta.
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
    pose.pdb_info = PDBInfo(pose)
    pose.pdb_info.name(seq[:8])

    return pose

# Monkey-patch read-only attributes
# Ex: Pose.size() becomes Pose.size
# Pose.size = 30 would return an error.
_read_attributes = ('const_data_cache', 'data', 'energies', 'atom_tree',
                    'membrane_info', 'num_chains', 'num_jump', 'observer_cache',
                    'reference_pose_set', 'reference_pose_set_cop', 'sequence',
                    'size', 'total_atoms', 'total_residue')

for attr in _read_attributes:
    setattr(Pose, attr, property(getattr(Pose, attr)))

# Create torsion tuple
# Ex: Pose.psis: tuple of psi.
# Pose.psis[5] returns Pose.psi(6)
# Pose.psis[5] = 180. uses Pose.set_psi(6, 180.)
_torsions_attributes = ('alpha', 'beta', 'delta', 'epsilon', 'gamma', 'mu',
                        'omega', 'theta', 'zeta', 'phi', 'psi')

for attr in _torsions_attributes:
    torsion_list_name = '{}s'.format(attr)
    getter = getattr(Pose, attr)
    setter = getattr(Pose, 'set_{}'.format(attr))
    setattr(Pose, torsion_list_name, torsion_list_property(getter, setter))


# Read/Set property
Pose.fold_tree = property(Pose.fold_tree, Pose.fold_tree)
Pose.pdb_info = property(Pose.pdb_info, Pose.pdb_info)
