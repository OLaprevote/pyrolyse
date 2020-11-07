# Dark magic, temporary if possible.
class _CallableProperty(property):
    def __init__(self, fcall=None, fget=None, fset=None, fdel=None, doc=None):
        self._call = fcall
        super().__init__(fget, fset, fdel, doc)

    def __call__(self, *args): return self._call(self, *args)


class AngleMaxDict(dict):
    """Dict where __setitem__ sets max_angle values from an object

    Parameters:
    -----------
    obj
        Object with an angle_max method.
    """
    def __init__(self, obj):
        self.obj = obj
        angles_dict = {sec: obj.get_angle_max(sec) for sec in 'HEL'}
        super().__init__(self, **angles_dict)

    def __setitem__(self, index, value):
        self.obj.angle_max(index, value)


def _get_angles_max(self):
    return(AngleMaxDict(self))


def _set_angles_max(self, angles):
    for sec in angles.keys():
        self.angle_max(sec, angles[sec])
