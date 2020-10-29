class TorsionList(list):
    def __init__(self, pose_instance, torsion_setter, *args):
        self.pose_instance = pose_instance
        self.torsion_setter = torsion_setter
        super().__init__(*args)

    def __setitem__(self, index, value):
        if isinstance(index, int):
            if index >= 0:
                self.torsion_setter(self.pose_instance, index+1, value)
            else:
                resid = self.pose_instance.size + index + 1
                self.torsion_setter(self.pose_instance, resid, value)

        # Slices are obtained when using semi-colons to get items, e.g. a[2:8:3]
        elif isinstance(index, slice):
            # Find concerned residue indexes
            residues = range(1, self.pose_instance.size+1)[index]
            for resid, new_psi in zip(residues, value):
                self.torsion_setter(self.pose_instance, resid, new_psi)

        else: raise TypeError(("TorsionList indices must be integers or"
                               " slices, not {}").format(type(index).__name__))


def torsion_list_property(getter, setter):
    def get_torsions(pose_instance):
        torsion_list = (getter(pose_instance, resid) for resid in range(1, pose_instance.size+1))
        return TorsionList(pose_instance, setter, torsion_list)

    def set_torsions(pose_instance, new_torsions):
        for resid in range(pose_instance.size):
            setter(pose_instance, resid+1, new_torsions[resid])

    return property(get_torsions, set_torsions)
