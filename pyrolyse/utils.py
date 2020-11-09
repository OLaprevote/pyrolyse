from io import StringIO
import logging
import os
from pathlib import Path

import pyrosetta as pyr
import pyrosetta.rosetta as ros

from pyrosetta.io import pose_from_file
from pyrosetta.toolbox.rcsb import pose_from_rcsb


__all__ = ['logger', 'get_pose', 'init', 'digest']


# WIP
logger = logging.getLogger('rosetta')
logger.setLevel(logging.WARNING)

# LOG catches logs, the idea is to have a variable you can call if you
# retrospectively want to see logs while stdout is disabled.
# Currently useless: catching anything only if rosetta logger level set
# to INFO. Also doesn't disable log to be printed to stdout.
#
# LOG = StringIO()
# stream_handler = logging.StreamHandler(LOG)
# stream_handler.setLevel(logging.INFO)


def pose_from_sequence(seq, res_type='fa_standard', auto_termini=True):
    pose = ros.core.pose.Pose()
    ros.core.pose.make_pose_from_sequence(pose, seq, res_type, auto_termini)

    try:
        total_residue = pose.total_residue()
        pose.pdb_info(ros.core.pose.PDBInfo(pose))
        pdb_info = pose.pdb_info()
    except TypeError:
        total_residue = pose.total_residue
        pose.pdb_info = ros.core.pose.PDBInfo(pose)
        pdb_info = pose.pdb_info

    for i in range(1, total_residue + 1):
        res = pose.residue(i)
        if not res.is_protein() or res.is_peptoid() or res.is_carbohydrate():
            continue

        pose.set_phi(i, 180)
        pose.set_psi(i, 180)
        pose.set_omega(i, 180)

    # Empty PDBInfo (rosetta.core.pose.PDBInfo()) is not correct here;
    # we have to reserve space for atoms....
    try:
        pdb_info.name(seq[:8])
    except TypeError:
        pdb_info.name = seq[:8]

    return pose


def get_pose(pose_input, kind=None):
    """Automatically get pose

    Guess what kind of pose fetcher to use based on the input.
    If it is a path object, will invoke pyrosetta.pose_from_file.
    Else will guess if it is a file or not, if not will either
    fetch in the rcsb or make from sequence.

    Parameters
    ----------
    pose_input: str or libpath.Path
        Either path to a file, RCSB ID or an amino-acid sequence.
    kind: {None, 'file', 'rcsb', 'sequence'}
        Which kind of input it is. Inferes which method to use from input
        if set to None.

    Output
    ------
    pose: pyrosetta.rosetta.core.pose.Pose or pyrolyse.pose.Pose
    """
    if kind is None:
        if Path(pose_input).exists():
            return pose_from_file(str(pose_input))

        elif any(c in str(pose_input) for c in ['.', os.sep]):
            raise FileNotFoundError(
                'No such file or directory: {}'.format(pose_input))

        elif len(pose_input) == 4 and any(c.isdigit() for c in pose_input):
            return pose_from_rcsb(pose_input)

        else:
            return pose_from_sequence(pose_input)

    else:
        if kind == 'file':
            return pose_from_file(str(pose_input))

        elif kind == 'rcsb':
            return pose_from_rcsb(pose_input)

        elif kind == 'sequence':
            return pose_from_sequence(pose_input)

        else: raise TypeError(("Argument kind can only be str 'file', "
                "'rcsb', 'sequence' or None type. Argument entered: {value}, "
                "{kind_type}").format(value=kind, kind_type=type(kind)))


def init(options='-ex1 -ex2aro', extra_options='', set_logging_handler=True,
         notebook=None, silent=False):
    """Init pyrosetta using pyrosetta logger at WARNING level

    To think about:
    - Would be better if Rosetta logger matched C++ log levels. Currently
    everything is put to INFO level.
    - The logs would go directly in io.StringIO object in case of Notebook.
    Only logs from warning or higher in stdout, the rest can be checked
    by calling pyrolyse.LOG. It can be chosen at first if StringIO used or
    outputted to a file, anyway StringIO can always be outputted afterward.
    Should find a way to "flush" StringIO object in case it becomes too big,
    or store it on disk to not have it in RAM. Tried to set-up a StringIO output
    through a StreamHandler but logger was also outputting to stdout.
    - In case of Notebooks, would probably be good to have a list-object,
    where one item appended after each use. This way it would be easy to call
    pyrolyse.LOG[-1] to simply see what happened during last cell used?
    """
    #TODO
    # Find way to totally redirect stdout

    pyr.init(options, extra_options, set_logging_handler, notebook, silent)


def digest():
    """Output a pandas DataFrame for any PyRosetta object
    """
    pass
