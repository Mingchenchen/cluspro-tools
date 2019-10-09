def interface(sele1, sele2, interface_radius=9.0):
    """
    Identifies the binding interface of a ligand-receptor complex
    :param sele1: (ProDy AtomGroup) with dimensions [[N, [xn, yn, zn]]
    :param sele2: (ProDy AtomGroup) with dimensions [[N, [xn, yn, zn]]
    :param interface_radius: float representing interface radius
    :return: coordinates of sele1 and sele2 interfaces
    """
    sele1_interface = sele1.select("all within {} of sele2".format(interface_radius), sele2=sele2)
    sele2_interface = sele2.select("all within {} of sele1".format(interface_radius), sele1=sele1)

    return sele1_interface, sele2_interface
