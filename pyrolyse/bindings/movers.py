# Dark magic, temporary if possible.
class _CallableProperty(property):
    """A callable property. property child.

    Parameters
    ----------
    fcall : function
        Function to be called by __call__
    fget : function
        __get__ method
    fset : function, optional
       __set__ method
    fdel : function, optional
        __del__ method
    doc : function, optional
        Documentation. If not imputed deduced from fget.

    See also
    --------
    pyrolyse.movers.simple.SmallMover.movemap
    pyrolyse.movers.simple.ShearMover.movemap
    """
    def __init__(self, fcall, fget, fset=None, fdel=None, doc=None):
        self._call = fcall
        super().__init__(fget, fset, fdel, doc)

    def __call__(self, *args): return self._call(self, *args)


class AngleMaxDict(dict):
    """Dict where __setitem__ sets max_angle values from an object

    Parameters
    ----------
    mover
        Object with an angle_max method, like SmallMover

    See also
    --------
    pyrolyse.movers.simple.SmallMover
    pyrolyse.movers.simple.ShearMover
    """
    def __init__(self, obj):
        self.obj = obj
        angles_dict = {sec: obj.get_angle_max(sec) for sec in 'HEL'}
        super().__init__(self, **angles_dict)

    def __setitem__(self, index, value):
        self.obj.angle_max(index, value)


def _get_angles_max(self):
    """Maximum angles of mover, in degrees

    Setting an item will modify the maximum angle.

    Examples
    --------
    >>> small_mv = lys.movers.simple.SmallMover()
    >>> small_mv.angles_max['E'] = 90.
    """
    return(AngleMaxDict(self))


def _set_angles_max(self, angles):
    """Setter of max_angles"""
    for sec in angles.keys():
        self.angle_max(sec, angles[sec])
