from io import StringIO
import logging

import pyrosetta as ros


__all__ = ['logger', 'LOG', 'get_pose', 'init', 'digest']


logger = logging.getLogger('rosetta')
log.setLevel(logging.WARNING)

# LOG catches LOG, the idea is to have a variable you can call if you
# retrospectively want to see logs while the stdout is disabled.
# Currently catching anything only if rosetta logger level set to INFO.
# Also doesn't disable log to be printed to stdout.
LOG = StringIO()

stream_handler = logging.StreamHandler(LOG)
strean_handler.setLevel(logging.INFO)


def get_pose(pose_path, new_pose=None):
    """Automatically get pose

    Guess what kind of pose fetcher to use based on the input.
    If it is a path object, will invoke pyrosetta.pose_from_file.
    Else will guess if it is a file or not, if not will either
    fetch in the rcsb or make from sequence.

    Attribute
    ---------
    pose_path: str or path
        Either path to a file, an amino-acid sequence or an RCSB ID.
    new_pose: bool or None
        Whether to return a pyrosetta (False) or a pyrolyse pose (True).
        Use pyrolyse.PYTHONIZE variable if None (False by default).

    Output
    ------
    pose: pyrosetta.rosetta.core.pose.Pose or pyrolyse.pose.Pose
    """
    pass


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

    ros.init(options, extra_options, set_logging_handler, notebook, silent)


def digest():
    """Output a pandas DataFrame for any PyRosetta object
    """
    pass
