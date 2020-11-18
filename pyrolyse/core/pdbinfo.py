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

PDBInfo.__doc__ = """Maintains pdb residue & atom information inside a Pose

Upon creation of new residue records, e.g. when calling the
constructors without 'init' or appending/prepending residues, the
chain letter for the new records will be set to a character, currently
'^', denoting "empty record".  This character may be looked up by
calling the static method PDBInfo::empty_record().

Class implementation is biased towards simplicity and fast lookup.
Residue/atom information are kept in vectors.  An internally maintained
PDBPoseMap provides mapping from pdb -> pose residue numbering. This
causes residue mutators to be a bit more expensive due to map updates,
but this is ok because they are typically called sparingly. Accessors
and mutators have overloaded method convention, while special mutators
use .set_* convention.

Modified by pyrolyse at pyrolyse.core.pdbinfo.

Parameters
----------
1. No parameter
2. n : int
3. pose : pyrosetta.rosetta.core.pose.Pose

4.
pose : pyrosetta.rosetta.core.pose.Pose
init : bool

5. arg0 : pyrosetta.rosetta.core.pose.PDBInfo

Examples
--------
>>> new_pdb_info = ros.core.pose.PDBInfo(pose)
>>> new_pdb_info.name = 'New name'
>>> pose.pdb_info = new_pdb_info   # If pyrolyse pose
"""
