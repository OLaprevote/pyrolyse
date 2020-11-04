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
# for attr in _torsions_attributes:
#     setattr(Pose, torsion_list_name, torsion_list_property(chi, set_chi))
