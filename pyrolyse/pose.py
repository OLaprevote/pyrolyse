"""Pyrosetta Pose and utilities monkey-patch

TODO

conformation and set_new_conformation: getter/setter?

Check residues, reslabels and scores.

append_residue: summarize all append_residue methods?
batch_get_xyz and batch_set_xyz: see how it works.
insert_residue: see if possible to summarize insert_residue_by_bond,
                insert_residue_by_jump and insert_residue_by_atoms.

chi: (chino: int, seqpos: int) a numpy array? -> Wouldn't work, not
     same dimension.
     Could create a list of list with tuples and tuples of slices
     comprehension.

aa: transform into a list?

pdb_rsd: pandas DataFrame-ish (takes a tuple (chain, resNo) as item)?
         like pdb_rsd['A', 10]? Possible to write pdb_rsd['A', 10:20],
         pdb_rsd['A'].
         Would pdb_rsd['A':'C', 10:20] work?

replace_residue     # To use as a residue setter for residues?

secstruct: pose.secstruct[24] returns secstruct of residue 24?
set_secstruct: then settable for range of residues
    However, pose.secstruct() returns the sequence of secondary
    structure while pose.secstruct(24) returns the secondary structure
    for said residue: maybe this behaviour shouldn't be broken.
"""

from pyrosetta.io import pose_from_sequence
from pyrosetta.rosetta.core.io.raw_data import ScoreMap
from pyrosetta.rosetta.core.simple_metrics import clear_sm_data
from pyrosetta.rosetta.core.pose import (Pose, clearPoseExtraScores,
                                         make_pose_from_sequence, PDBInfo)
from pyrosetta.bindings.pose import (PoseResidueAccessor,
                            PoseResidueLabelAccessor, PoseScoreAccessor)

from .pythonizer.pose import TorsionList, torsion_list_property


__all__ = ['Pose', 'pose_from_sequence']


def pose_from_sequence(seq, res_type="fa_standard", auto_termini=True):
    """Returns a pose from single-letter protein sequence.

    Modified by pyrolyse from pyrosetta.
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


# Patch data descriptors methods calling monkey-patched attribute.
def _len_residues(self):
    return self.pose.size


PoseResidueAccessor.__len__ = _len_residues


def _len_reslabels(self):
    return self.pose.pdb_info.nres()


def _getitem_reslabels(self, key):
    """1-based index and slice over residue labels."""
    if isinstance(key, slice):
        return (self[i] for i in range(*slice_1base_indicies(key, len(self))))
    else:
        if key == 0:
            raise IndexError("1 base indexing does not support 0 index")
        if key < 0:
            key = len(self) + 1 + key
        return ResidueLabelAccessor(self.pose.pdb_info, key)


PoseResidueLabelAccessor.__len__ = _len_reslabels
PoseResidueLabelAccessor.__getitem__ = _getitem_reslabels


def _get_energies_scores(self):
    import types

    return types.MappingProxyType(self.pose.energies.active_total_energies())


def _get_all_scores(self):
    import types

    return types.MappingProxyType(
        dict(
            list(ScoreMap.get_arbitrary_string_data_from_pose(self.pose).items())
            + list(ScoreMap.get_arbitrary_score_data_from_pose(self.pose).items())
            + list(self.pose.energies.active_total_energies().items())
        )
    )


def _clear_scores(self):
    """ Clear pose energies, extra scores, and SimpleMetric data"""
    self.pose.energies.clear()
    clearPoseExtraScores(self.pose)
    clear_sm_data(self.pose)


PoseScoreAccessor.energies = property(_get_energies_scores)
PoseScoreAccessor.all = property(_get_all_scores)
PoseScoreAccessor.clear = _clear_scores
