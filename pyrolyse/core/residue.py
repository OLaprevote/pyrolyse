"""Pyrosetta monkey-patched Residues

TODO
----
chi -> returns chi torsion angles
set_chi -> item setter
set_all_chi as setter.
"""
from pyrosetta.rosetta.core.conformation import Residue


__all__ = [Residue]

# Monkey-patch methods to read-only attributes
# Ex: Residue.n_nus() becomes Residue.n_nus
# Pose.size = 30 would return an error.
_read_attributes = ('aa', 'accpt_pos', 'accpt_pos_sc', 'actcoord',
    'actcoord_atoms', 'all_bb_atoms' , 'atoms', 'atoms_with_orb_index',
    'backbone_aa', 'carbohydrate_info', 'data',
    'last_backbone_atom', 'lower_connect', 'lower_connect_atom',
    'mainchain_atoms', 'n_hbond_acceptors', 'n_hbond_donors',
    'n_mainchain_atoms', 'n_non_polymeric_residue_connections', 'n_nus',
    'n_orbitals', 'n_polymeric_residue_connections',
    'n_possible_residue_connections', 'n_virtual_atoms', 'na_analogue', 'name',
    'name1', 'name3', 'natoms', 'nbr_atom', 'nbr_radius', 'nchi',
    'nheavyatoms', 'nonconst_data', 'path_distances','type',
    )
for attr in _read_attributes:
    setattr(Residue, attr, property(getattr(Residue, attr)))

# Read/Set property
_read_write_attributes = ('mainchain_torsions', 'seqpos',)
for attr in _read_write_attributes:
    getter_setter = getattr(Residue, attr)
    setattr(Residue, attr, property(getter_setter, getter_setter))

# TODO
# Create chi torsion tuple

Residue.__doc__ = """Instance Residue class, used for placed residues and rotamers

This class is designed to be lightweight.
It holds a const-reference ("rsd_type_") to a ResidueType object for
access to information common to all instances of a single type, e.g.,
Alanine or Thymine.  Residue stores any data unique to a placed residue
or rotamer, currently:

- a vector1 of Atoms, which holds the positions (and also the atom-types
for fast access during scoring);

- the sequence position and chain number, both integers
(See the documentation of Pose::num_chains() for details about chain
numbers, chain letters and jumps.)

- the backbone, side-chain, and internal ring (if applicable) torsion
angles (of course backbone torsions are not unique to a rotamer, and
the chi angles are derivable from the coordinates, but storing them in
the residue is convenient for scoring purposes).

- the coordinates of an interaction center or centroid, used e.g., in the
knowledge-based full-atom pair term ("actcoord_").  Maybe this will also
hold the centroid position for centroid-mode scoring??

Modified by pyrolyse in pyrolyse.core.residue

Parameters
----------
1.
rsd_type_in: pyrosetta.rosetta.core.chemical.ResidueType
dummy_arg: bool

2.
rsd_type_in: pyrosetta.rosetta.core.chemical.ResidueType
dummy_arg: bool

3.
rsd_type_in: pyrosetta.rosetta.core.chemical.ResidueType
current_rsd: pyrosetta.rosetta.core.conformation.Residue
conformation: pyrosetta.rosetta.core.conformation.Conformation

4.
rsd_type_in: pyrosetta.rosetta.core.chemical.ResidueType
current_rsd: pyrosetta.rosetta.core.conformation.Residue
conformation: pyrosetta.rosetta.core.conformation.Conformation
preserve_c_beta: bool

5.
rsd_type_in: pyrosetta.rosetta.core.chemical.ResidueType
current_rsd: pyrosetta.rosetta.core.conformation.Residue
conformation: pyrosetta.rosetta.core.conformation.Conformation
preserve_c_beta: bool
allow_alternate_backbone_matching: bool

6. arg0: pyrosetta.rosetta.core.conformation.Residue

Attributes
----------
aa
natoms
type

Examples
--------
>>> res = pose.residue(1)
>>> print(res.natoms)
>>> print(res.atom_index('CA'))
"""
