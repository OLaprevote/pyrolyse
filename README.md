# Pyrolyse, a wrapper for PyRosetta

Pyrolyse is a python wrapper around PyRosetta, which is already a python
wrapper around Rosetta, written in C++.
So why wrapping? Because currently, PyRosetta feels like a C++ interface.
Docstrings aren't this informative if you don't know how they are written
in C++, there is some boilerplate code, the organization doesn't seem
quite logical either, and a lot of objects are Rosetta types while they
feel like a simple dictionary, a numpy array or a pandas dataframe is
all you'd need to feel right at ease.

Pyrolyse NEEDS PyRosetta to be installed.

Hence Pyrolyse is intended to be a digest of Pyrosetta, but they can
absolutely work hand in hand. They are actually intended to work hand
in hand: this way you can take best of both worlds!

Here are a few things I'd want to achieve:

1. The standard output of most pyrosetta functions (like `init`,
`pose_from_xxx` or `get_function_blabla` always give a standard output.
I personnally don't care about it, so I would like to either store it
in a global variable pyrolyse.LOG. Still outputing errors, though.

2. Run `pyrosetta.init()` directly during import. Still possible to rerun
it afterward if some special flags are needed, or put said flags in
`./rosetta/flags`, but if you import pyrosetta anyways you will `init()`
so why the heck do I have to write it every time?
Also, if first point seems too much of a hassle, `init(silent=True` and with
options to show only warnings level logs would be nice.
`pyrolyse.init` function should have a `log_to_glob: bool` argument.

3. Have a `pyrolyse.digest` function. For PyRosetta objects like
`Energies`, `Pose`, `Residue` or `Pose.residues` would automatically
output a `pandas DataFrame` with the important infos. For other objects
would output a dictionary, a list or a `numpy array`. Of course there
would be several functions like `pyrolyse.list`, `pyrolyse.dict`,
`pyrolyse.dataframe`, etc and `pyrolyse.digest` would apply given
different types. Should test before if just a dumb `list` doesn't work, though.

4. Transform PyRosetta objects to a more pythonic form.
For example, `Residue` objects have many methods which should simply be
attributes, like `residue.name()` which will never take an argument,
so it could actually be `residue.name`. using `@property` decorator.
In the same idea, in the case of `residue.chi()` could be `@property` and
transformed in a `numpy.array`, this way `residue.chi[0] = 5.` could be
done, instead of the unpythonic `residue.set_chi(1, 5.)` (which
would still work). There could also be a `@chi.setter` which would verify
that the length of the array would be good if one suddenly wanted to do
`residue.chi = [180.5, 25.0, 10.]`. This behavior would allow easy
ways to rotate everything from different angles at once, like
`residue.chi[:2] += [25., 10.]` if I wanted to rotate the first chi angle
of 25° and the second one of 10°. While this wouldn't probably be used in
a serious PyRosetta script, the behavior would still feel more natural.

5. Make the variable locations actually make sense, for example with
objects for xml import or Movers.
Why write
```
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
min_mover = MinMover()
```
When I could write
```
from pyrolyse import movers
min_mover = movers.Min()
```

6. Automatize
Take PyRosetta's `pose_from_xxx`. There should be a dumb `get_pose` which
simply identifies wich `pose_from_xxx` object is needed given the input.
(If it's a path it WILL be a file, if a string first check if it has file
extension, if not check if it is only 4 letters (-> probably rcsb, check if
it has numbers or has only lower case). Alternatively there should be a
`from` argument (default `None`) where one could specify the type of input.
This is true for Movers as well: why write
```
from pyrosetta.teaching import MinMover
sfxn = foo1()
movemap = foo2()
min_mover = MinMover()
MinMover.set_score_function(sfxn)
MinMover.set_movemap(movemap)
```
when I could write
```
from pyrolyse import movers
sfxn = foo1()
movemap = foo2()
min_mover = movers.Min(score_function=sfxn, movemap=movemap)
```
There should also be a function `apply_movers` wich simply takes
a pose and a list/tuple of movers, then apply them sequencially to said pose.
In the same tone, the `SequenceMover` should also be able to take a list of movers.

7. Use XML as functions
RosettaScript allows variables, yet it didn't seem set in PyRosetta.
There should be a simple `xml_func` taking a xml with variables in it
and creating a function from it. This would allow the use of much MUCH more
interesting use of PyRosetta: the xml used for peptide generation in 2016 article
from Chevalier & al. could be distributed as a simple PyRosetta function.
