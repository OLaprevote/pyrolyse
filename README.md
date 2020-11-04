# Pyrolyse, syntactic sugar for PyRosetta

Pyrosetta is a python API of Rosetta, a software intended for protein
design.

**Pyrolyse NEEDS PyRosetta to be installed.**

Parts of PyRosetta source code are copied and modified here.
License for PyRosetta can be found here:
<https://els2.comotion.uw.edu/print/licencing-option/19>

Pyrolyse mostly performs monkey-patching of existing classes in order to make
PyRosetta feel more pythonic. There are also some dynamic definitions of new
class attributes, and some useful new functions are created as well, 

The initiative in itself is probably shallow, and was simply started because
some PyRosetta methods didn't feel "palatable" to me, which gives you an idea of
how deep the stick must be up my [Argininosuccinate Synthase][1].

As PyRosetta and Rosetta are ever-growing softwares, not everything can be
patched manually so that it has a "pythonic" vibe.
Which is why having an independant module doing it may be a good idea:
while not needed to run PyRosetta, it can be imported over it to smoothen
its use where the smoothing has been developed. If it is not yet patched
or buggy, one can always use the unpatched pyrosetta function.

For example, in pyrosetta, to get psi of the first residue then set it to zero
one would write:
```python
print(pose.psi(1))
pose.set_psi(1, 0.)
```

While when using pyrolyse over pyrosetta, one could also write:
```python
print(pose.psis[0])
pose.psis[0] = 0.
```

Which also permits to redefine several psis in a row with slices (although directly
modifying torsion angles may not be the first goal when using pyrosetta).
Other simple things consist in transforming most methods serving as read-only
attributes into, well, read-only attributes, e.g. `pose.size()` becomes
`pose.size`. This can be useful during autocompletion in notebooks or python
interpreter to quickly grasp the usage of each proposition.

[1]:https://www.uniprot.org/uniprot/P00966

# How to install

Given licensing problems as some functions are copied (although
modified) from PyRosetta, pyrolyse was not made available on PyPI for
easy `pip install`. Hence the best way to get it would be to get this
repository on your computer, go inside, then run `pip install -e .`.
Note that the environment in which it is installed must contain PyRosetta
package.

# How to run

There are multiple ways to use pyrolyse:
 - Either keep importing everything from PyRosetta, then import `pyrolyse.all`
   to monkey-patch every objects and import new defined functions, like:

```Python
import pyrosetta as pyr
import pyrosetta.rosetta as ros
import pyrolyse.all as lys


lys.init() # Same as pyr.init excep set_logging_handler set to True.
pose_lys = lys.get_pose('5WRG')  # Gets structure from rcsb
print(pose_lys.size, pose_lys.residue(1).n_chis) 

# Former pyrosetta functions will also output monkey-patched classes
pose_pyr = pyr.pose_from_pdb('5WRG.pdb')
print(pose_pyr.size)

try:
    pose_error = pyr.pose_from_sequence('A'*12)
except:
    print('Sadly, some functions from pyrosetta crash')

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

pose = pyr.pose_from_pdb('
```

  - Or to import some useful functions without changing pyrosetta behavior:

```python
from pyrolyse.utils import get_pose

pose = get_pose('LYSE')  # Recognize it as a sequence
```

 - As one of the end-term goal would be to cover most objects from PyRosetta,
   it should also be possible to work while importing mostly pyrolyse objects
   which should feature a less intricate way to access PyRosetta objects, at
   least for protocols creation:

```python
import pyrolyse.all as lys
from pyrolyse.movers import simple
from pyrosetta import MoveMap

lys.init()
lys.logger.setLevel('INFO')

pose = lys.get_pose('A'*12)
mmap = MoveMap()
mmap.set_bb(True)

small_mover = simple.SmallMover(movemap, 1., 1)
small_mover(pose)   # Same than small_mover.apply
```

More examples are showcased in the notebook directory while following PyRosetta notebooks
tutorials.

# Work In Progress

## Logging
Make python logging levels fit with Rosetta C++ ones. Currently everything
from info to warnings is outputed to logging.INFO level.

Setup pyrolyse to automatically stdout only warnings or higher: anyway the
logger can be set to `logging.INFO` with simple
`pyrolyse.logger.setLevel('INFO')`.

Just in case, store in a module variable `pyrolyse.LOG` which registers
some of the last logs even at info level. This way if a particularly long
process was executed and there is a problem not showing in warnings, one
can get the last info logs from pyrolyse.LOG.

Enable the above behavior only if python is interactive (console or notebook),
NOT in script: in this case everything at INFO or less (depending on flags
in init) should be shown.

# TODO

## Digest
Have a `pyrolyse.digest` function. For PyRosetta objects like
`Energies`, `Pose`, `Residue` or `Pose.residues` would automatically
output a `pandas DataFrame` with the important infos. For other objects
would output a dictionary, a list or a `numpy array`. There
could be several functions like `pyrolyse.list`, `pyrolyse.dict`,
`pyrolyse.dataframe`, etc and `pyrolyse.digest` would apply them given
different types.

## Make movers and movemap setup easier
Be able to write:
```python
import pyrolyse.all as lys

sfxn = foo1()
movemap = lys.MoveMap(bb=True)
min_mover = lys.movers.simple.SmallMover(score_function=sfxn, movemap=movemap)
```

instead of:
```python
from pyrosetta.teaching import SmallMover

sfxn = 
mmap = MoveMap()
mmap.set_bb(True)
small_mover = SmallMover(mmap, 1., 1)
MinMover.set_score_function(sfxn)
```

Add a function `apply_movers` wich takes a pose and a generator of movers,
then apply them sequencially to said pose.
Although it is possible to do it currently in one line with:
```python
map(lambda mover: getattr(mover, apply)(pose), list_movers)
```

Make `SequenceMover` able to take in a list of movers.

## Use XML as functions
RosettaScript allows variables in xml, yet it didn't seem set in PyRosetta.
(Though I haven't investigated yet.)
Add a simple `xml_func` taking a xml with variables in it
and creating a function from it, like:
```python
new_pyrosetta_protocol = pyrolyse.xml_func(xml_string)
new_pyrosetta_protocol(protein_to_dock, number_of_protein_to_design)
```

Case example: xml used for peptide generation in 2016 article
from Chevalier & al. could be distributed as a simple PyRosetta function.

