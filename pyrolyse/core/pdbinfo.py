"""Pyrosetta monkey-patched PDBInfo

TODO
----
chain: settable list
bfactor: settable list
set_chains set_crystinfo
"""
from pyrosetta.rosetta.core.pose import PDBInfo


__all__ = [PDBInfo]

# Monkey-patch methods to read-only attributes
# Ex: Residue.n_nus() becomes Residue.n_nus
# Pose.size = 30 would return an error.
_read_attributes = ('nres', 'num_chains', 'short_desc',)
for attr in _read_attributes:
    setattr(PDBInfo, attr, property(getattr(PDBInfo, attr)))

# Read/Set property
_read_write_attributes = ('header_information', 'modeltag', 'name', 'obsolete',
    'remarks',
    )
for attr in _read_write_attributes:
    getter_setter = getattr(PDBInfo, attr)
    setattr(PDBInfo, attr, property(getter_setter, getter_setter))

_get_set_pairs = (('crystinfo', 'set_crystinfo'),)
for get_attr, set_attr in _get_set_pairs:
    getter = getattr(PDBInfo, get_attr)
    setter = getattr(PDBInfo, set_attr)
    setattr(PDBInfo, get_attr, property(getter, setter))
