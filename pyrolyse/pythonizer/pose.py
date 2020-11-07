class TorsionList(list):
    def __init__(self, pose, torsion_setter, *args, **kwargs):
        self.pose = pose
        self.torsion_setter = torsion_setter
        super().__init__(*args, **kwargs)

    def __setitem__(self, index, value):
        if isinstance(index, int):
            if index >= 0:
                self.torsion_setter(self.pose, index+1, value)
            else:
                resid = self.pose.size + index + 1
                self.torsion_setter(self.pose, resid, value)

        # Slices are obtained when using semi-colons to get items, e.g. a[2:8:3]
        elif isinstance(index, slice):
            # Find concerned residue indexes
            residues = range(1, self.pose.size+1)[index]
            for resid, new_psi in zip(residues, value):
                self.torsion_setter(self.pose, resid, new_psi)

        else: raise TypeError(("TorsionList indices must be integers or"
                               " slices, not {}").format(type(index).__name__))


def torsion_list_property(getter, setter):
    # TODO: Change this for DNA torsion lists
    def get_torsions(pose):
        torsion_list = (getter(pose, resid) for resid in range(1, pose.total_residue+1)
                        if any((pose.residue(resid).is_protein(),
                                pose.residue(resid).is_carbohydrate(),
                                pose.residue(resid).is_peptoid()
                                ))
                        )
        return TorsionList(pose, setter, torsion_list)

    def set_torsions(pose, new_torsions):
        for resid in range(pose.size):
            setter(pose, resid+1, new_torsions[resid])

    return property(get_torsions, set_torsions)
