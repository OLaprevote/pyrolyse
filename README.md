# Pyrolyse, syntactic sugar for PyRosetta

PyRosetta is a python API of Rosetta, a software intended for protein
design.

**Pyrolyse NEEDS PyRosetta to be installed.**

This package contains modified PyRosetta source code.
Academic License for PyRosetta can be read here:
<https://els2.comotion.uw.edu/print/licencing-option/19>

Pyrolyse mostly performs monkey-patching of existing classes in order to make
PyRosetta feel more pythonic. There are also some dynamic definitions of new
class attributes, and some useful new functions are created as well, 

For example, in PyRosetta, to get psi of the first residue then set it to zero
one would write:
```python
print(pose.psi(1))
pose.set_psi(1, 0.)
```

While when using pyrolyse over PyRosetta, one could also write:
```python
print(pose.psis[0])
pose.psis[0] = 0.
```

Which also permits to redefine several psis in a row with slices (although directly
modifying torsion angles may not be the first goal when using PyRosetta).
Other simple things consist in transforming most methods serving as read-only
attributes into, well, read-only attributes, e.g. `pose.size()` becomes
`pose.size`. This can be useful during autocompletion in notebooks or python
interpreter to quickly grasp the usage of each proposition.

As PyRosetta and Rosetta are huge ever-growing softwares, not everything can be
patched manually right away so that it gets a "pythonic" vibe.
Which is why having an independant module doing it may be a good idea:
while not needed to run PyRosetta, it can be imported over it to smoothen
its use, when it has been developed. If it is not yet patched
or too buggy, one can always use the unpatched PyRosetta function.

[1]:https://www.uniprot.org/uniprot/P00966

# How to install

Pyrolyse has not yet been made available on PyPI for
easy `pip install`. Hence the best way to get it would be to get this
repository on your computer using `git clone https://github.com/OLaprevote/pyrolyse.git`.
Go inside, then run `pip install -e .`.
Note that the environment in which it is installed must contain PyRosetta
package, and setuptools to do `pip`.  These are currently the only requirements.

# How to run

There are multiple ways to use pyrolyse:
 - Keep importing everything from PyRosetta, then import `pyrolyse.all`
   to monkey-patch every objects and import new defined functions, like:

```Python
import pyrosetta as pyr
import pyrosetta.rosetta as ros
import pyrolyse.all as lys


lys.init() # Same as pyr.init excep set_logging_handler set to True.
pose_lys = lys.get_pose('5WRG')  # Swiss-army function to get pose
print(pose_lys.size, pose_lys.residue(1).natoms) 

# Former PyRosetta functions will also output monkey-patched classes
pose_pyr = pyr.pose_from_pdb('5WRG.pdb')
print(pose_pyr.size)

try:
    pose_error = pyr.pose_from_sequence('A'*12)
except:
    print('Sadly, some functions from PyRosetta crash')

# Although they should have an equivalent in pyrolyse
# For example either:
pose_seq = lys.get_pose('A'*12)
# Or:
pose_seq = lys.pose_from_sequence('A'*12)
```

 - It is also possible to cherry-pick which classes you want to monkey-patch:

```Python
import pyrosetta as pyr
import pyrolyse.pose

pyr.init()

pose = pyr.pose_from_pdb('5WRG.pdb')
print(pose.size)
print(pose.residue(1).natoms())  # Residue class not monkey-patched
```

  - Or to import some useful functions without changing PyRosetta behavior:

```python
from pyrolyse.utils import get_pose

pose = get_pose('LYSE')  # Recognize it as a sequence
```

More examples are showcased in PyRosetta notebook directory. 

## Todo

- Add missing docstrings.
- Add a line like "Modified by pyrolyse in pyrolise.movers.simple" in every docstring
  so that modification is easily seen and source code is easily found.
- Write a setup.py with more than one line and add a version following [semantic versioning](https://semver.org/).
- Then add a version tag on git. Or go clever parmed way and use the Versioneer.
- Patch Mover base class: hopefully its changes repercut on other movers.
- Probably not, though, but argument list can then be imported and concatenated
  for each movers.
- Finish to process movers which were partially but not fully monkey-patched.
- Modify Movemap attributes.
- Enhance logging.
- Make Digest function.
- Stop writing explanatory notebooks and convert the rest of PyRosetta.
- Check if all python functions of PyRosetta work with current changes (most probably not).
- Not important and maybe not a good idea: add torsion angles to residues, like `pose.residue(1).psi`.
