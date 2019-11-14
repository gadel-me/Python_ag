
import re
import math
import md_stars as mds

"""
Helper functions for ag_prmtop module.
"""

"""
#### Functions only for use within 'read_prmtop' ####
"""


def chunk_list(lst, n):
    """
    Divide a list (listname) into chunks with length n
    """
    return [lst[i:i + n] for i in range(0, len(lst), n)]


def parse_section(opened_prmtop, geps, chunksize=None, itype=None):
    """
    Helper function for parsing an amber prmtop-file.
    geps:   given entries per section
    mepl:   max possible entries per line
    cpe:    characters per entry
    mepl:   max entries per line (entries given by format line in prmtop)
    lps:    lines per section (e.g. AMBER_ATOM_TYPE)
    """
    line      = next(opened_prmtop)  # line with formatting info
    mepl, cpe = [int(i) for i in re.findall(r'\d+', line)][:2]
    lps       = int(math.ceil(geps/mepl))
    entries   = []

    if "20a4" in line:
        # following lines do not necessarily have to have spaces
        for iline in range(lps):
            line = opened_prmtop.next().rstrip("\n")  # remove next line
            entries.extend([line[i:i+cpe].strip() for i in
                            range(0, len(line), cpe)])
    else:
        # following lines definitely have spaces and consist of floats or ints
        for iline in range(lps):
            line = next(opened_prmtop)
            if itype == "int":
                entries.extend([int(i) for i in line.split()])
            else:
                entries.extend([float(i) for i in line.split()])

    if chunksize:
        return chunk_list(entries, chunksize)
    else:
        return entries


def unmask_imp(section_dih_inc_h, section_dih_non_h):
    """
    Amber does distinguish between imp and dih insofar that only a dihedral-
    style is defined (but none for improper-angles). But the prmtop different-
    iates them nonetheless marking the last atom of each dihedral with a sign.
    -> unsigned dihedrals are dihedrals, signed dihedrals are impropers.
    To not loose this information, a helper-dictionary is created so that imp-
    propers can still be distinguished when processing the dihedral styles.
    """
    dih_imp_dict = {}
    both = section_dih_inc_h + section_dih_non_h

    for dih in both:
        if dih[3] < 0:  # last atom of dihedral
            dih_imp_dict[dih[4]] = "imp"
        else:
            dih_imp_dict[dih[4]] = "dih"

    return dih_imp_dict


def check_cvff_compatibility(prm_d, prm_n):
    """
    Check if lammps convention for impropers is fulfilled:
        cos(prm_d) = -1.0 or 1.0
        prm_n in = 1 or 2 or 3 or 4 or 6
    Sources:    Reading: http://lammps.sandia.gov/doc/improper_cvff.html
    """
    cos_prm_d = math.cos(prm_d)

    # check if difference between cos(prm_d) and -1/+1 is not less than 1e-10
    if not(cos_prm_d+1 < 1e-10) or not(cos_prm_d-1 < 1e-10):
        print("***Warning: Improper not compatible with lammps improper_style 'cvff'!")
        print("            cos(prm_d) must be -1 or 1, not {:.4f}!".format(cos_prm_d))
    if prm_n not in (1, 2, 3, 4, 6):
        print("***Warning: Improper not compatible with lammps improper_style 'cvff'!")
        print("            prm_n should only be (1, 2, 3, 4, 6)!")


def true_atm_id(prmtop_atm_id):
    """
    Atom indexes are actually indexes into a coordinate array, so the
    actual atom index A is calculated from the coordinate array index N by
    A = N/3 + 1. (N is the value in the topology file)
    """
    # only positive integers
    prmtop_atm_id = abs(prmtop_atm_id)
    return int(prmtop_atm_id/3 + 1)  # true bond index


def relocate_parsed(sect_bad_iwh, container, key, dict_key_old_new, dict_atm_id_old_new):
    """
    Sort the parsed data from the sections:
        > bonds (inc/without hydrogen) or
        > angles (inc/without hydrogen) or
        > dihedrals (inc/without hydrogen)
    to their respective containers (bonds, angles, dihedrals; see class Universe
    for more information).
    'bad' = b(onds), a(ngles), d(ihedrals)
    sect_bad_iwh:  list; section inclusive hydrogens or section without hydrogens
    key:           "bonds", "angles", "dihedrals"
    """
    #lc = len(container)  # number of elements in container

    # bonds with hydrogens
    for ientry in sect_bad_iwh:
        # transform amber-numbered atom-ids to their original ids
        ientry_0 = true_atm_id(ientry[0])
        ientry_0 = dict_atm_id_old_new[ientry_0]

        ientry_1 = true_atm_id(ientry[1])
        ientry_1 = dict_atm_id_old_new[ientry_1]

        if len(ientry) > 3:  # angles or dihedrals
            ientry_2 = true_atm_id(ientry[2])
            ientry_2 = dict_atm_id_old_new[ientry_2]

        if len(ientry) > 4:  # dihedrals
            ientry_3 = true_atm_id(ientry[3])
            ientry_3 = dict_atm_id_old_new[ientry_3]

        if key == "bonds":
            container.append(mds.Bond(bnd_id=len(container),
                                      atm_id1=ientry_0,
                                      atm_id2=ientry_1,
                                      bnd_key=dict_key_old_new[ientry[2]]))
        elif key == "angles":
            container.append(mds.Angle(ang_id=len(container),
                                       atm_id1=ientry_0,
                                       atm_id2=ientry_1,
                                       atm_id3=ientry_2,
                                       ang_key=dict_key_old_new[ientry[3]]))
        elif key == "dihedrals":
            # only append non-impropers, which are dihedrals
            if ientry[3] > 0:
                container.append(mds.Dihedral(dih_id=len(container),
                                              atm_id1=ientry_0,
                                              atm_id2=ientry_1,
                                              atm_id3=ientry_2,
                                              atm_id4=ientry_3,
                                              dih_key=dict_key_old_new[ientry[4]]))
        elif key == "impropers":
            #print(ientry_0, ientry_1, ientry_2, ientry_3, ientry[4])
            # if third int of each subcontainer is < 0, amber assumes it is an improper
            if ientry[3] < 0:
                container.append(mds.Improper(imp_id=len(container),
                                              atm_id1=ientry_0,
                                              atm_id2=ientry_1,
                                              atm_id3=ientry_2,
                                              atm_id4=ientry_3,
                                              imp_key=dict_key_old_new[ientry[4]]))
        else:
            raise RuntimeError("key must be 'bonds', 'angles', 'dihedrals' or 'impropers'!")

    return container


def fix_keys(container, dict_w_ptrs, key):
    """
    Function 'relocate_parsed' makes some ordering, but cannot know which
    dihedral-/improper-keys are up-to-date (dih-/imp-types are renumbered while
    read consecutively). This function helps closing this gap. Keep in mind
    this function here is ONLY to HELP ag_prmtop-read_prmtop-function, since
    it is highly specific!
    """
    if key == "dihedrals":
        for dih in container:
            old_key = dih.dih_key
            new_key = dict_w_ptrs[old_key]
            dih.dih_key = new_key
    elif key == "impropers":
        for imp in container:
            old_key = imp.imp_key
            new_key = dict_w_ptrs[old_key]
            imp.imp_key = new_key
    else:
        raise KeyError("***'key' must be 'dihedrals' or 'impropers'!")


def entry_not_read(entry_name):
    """
    Helper function to print which entries are skipped.
    """
    print("***Prmtop-Info: Entry {} currently not read. Skipping.".format(entry_name))


def get_AB_ids(ntypes, section_nonbonded_parm_index):
    """
    Get ii-indices for acoef/bcoef-section.
    ii = Lennard-Jones-interactions for same atom type.
    Formula for acoef/bcoef-index: cur_ii_nonbon_param_id = ntypes*(cid-1)+cid
    """
    ii_nonbon_param_ids = []

    # 1. get lj-ii-indices for section 'nonbonded_parm_index'
    for i_id in range(ntypes):
        i_id += 1  # counting in prmtop starts with 1 (Fortran code)
        cur_ii_nonbon_param_id = ntypes*(i_id-1) + i_id
        # since we are not using Fortran, counting starts on 0
        # -> substract 1 from index at the end
        ii_nonbon_param_ids.append(cur_ii_nonbon_param_id-1)

    # 2. get lj-ii-indices for sections 'a_coef' and 'b_coef' from section
    #    nonbonded_parm_index
    ii_abcoef_ids = []
    for i in ii_nonbon_param_ids:
        # since we are not using Fortran, counting starts on 0 (-> subtract 1)
        cur_abcoef_id = section_nonbonded_parm_index[i] - 1
        ii_abcoef_ids.append(cur_abcoef_id)

    return ii_abcoef_ids


def sig_eps_from_AB(acoef, bcoef):
    """
    AB form: A = 4 * eps * sigma**12,
             B = 4 * eps * sigma**6
             sigma = (A/B)**(1/6)
             eps   = B**2/(4*A)
    Sources:    https://en.wikipedia.org/wiki/Lennard-Jones_potential#AB_form
    """
    try:
        sigma = (acoef/bcoef)**(1/6)
    except(ZeroDivisionError):
        sigma = 0.00000000E+00

    try:
        epsilon = bcoef**2/(4*acoef)
    except(ZeroDivisionError):
        epsilon = 0.00000000E+00

    return (sigma, epsilon)
