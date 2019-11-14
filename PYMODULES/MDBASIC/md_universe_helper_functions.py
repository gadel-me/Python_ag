
import ag_vectalg as agv
import ag_cryst as agc


def get_lattice(box):
        # get lattice box vectors for box
    if box.boxtype == "cartesian":
        box.cart2lat()
    elif box.boxtype == "lammps":
        box.lmp2lat()
    else:
        pass


def _crds_cart2fract(a, b, c, alpha, beta, gamma, cart_coords):
    """
    Convert cartesian coordinates to fractional coordinates.
    """
    fract_coords = []

    M_c2f = agc.M_cart2fract(a, b, c, alpha, beta, gamma)
    for coords in cart_coords:
        cur_fcoords = agv.mv_mult(M_c2f, coords)
        fract_coords.append(cur_fcoords)

    return fract_coords


def _del_nones(ilist):
    return [i for i in ilist if i is not None]


def merge_entries(dict1, dict2):
    """
    Merge two dictionaries and append not found matches. Helper function
    for _merge_dicts function.
    """
    num_entries = len(dict1)
    ptr_dict    = {}

    for key1 in dict1:
        for key2 in dict2:
            same = True

            # check if other attributes (where present) are different
            try:
                if dict1[key1].weigh != dict2[key2].weigh:
                    same = False
            except AttributeError:
                pass

            try:
                if dict1[key1].sitnam != dict2[key2].sitnam:
                    same = False
            except AttributeError:
                pass

            try:
                if dict1[key1].epsilon != dict2[key2].epsilon:
                    same = False
            except AttributeError:
                pass

            try:
                if dict1[key1].prm1 != dict2[key2].prm1:
                    same = False
            except AttributeError:
                pass

            try:
                if dict1[key1].prm2 != dict2[key2].prm2:
                    same = False
            except AttributeError:
                pass

            try:
                if dict1[key1].prm3 != dict2[key2].prm3:
                    same = False
            except AttributeError:
                pass

            try:
                if dict1[key1].prm4 != dict2[key2].prm4:
                    same = False
            except AttributeError:
                pass

            try:
                if dict1[key1].prm_k != dict2[key2].prm_k:
                    same = False
            except AttributeError:
                pass

            try:
                if dict1[key1].prm_n != dict2[key2].prm_n:
                    same = False
            except AttributeError:
                pass

            try:
                if dict1[key1].prm_d != dict2[key2].prm_d:
                    same = False
            except AttributeError:
                pass

            try:
                if dict1[key1].weigh_factor != dict2[key2].weigh_factor:
                    same = False
            except AttributeError:
                pass

            try:
                if dict1[key1].elec_inter_1_4_scale != dict2[key2].elec_inter_1_4_scale:
                    same = False
            except AttributeError:
                pass

            try:
                if dict1[key1].vdw_inter_1_4_scale != dict2[key2].vdw_inter_1_4_scale:
                    same = False
            except AttributeError:
                pass

            try:
                if dict1[key1].energy_unit != dict2[key2].energy_unit:
                    print("***Warning: Different energy units in entry. Please revise!")
            except AttributeError:
                pass

            if same is True:
                ptr_dict[key2] = key1
                dict2[key2].found = True
                break

    for key2 in dict2:
        try:
            if dict2[key2].found is True:
                pass
        except AttributeError:
            num_entries += 1
            dict1[num_entries] = dict2[key2]
            ptr_dict[key2] = num_entries

    return ptr_dict


# Key-functions for sorted() ---------------------------------------------------
# Sources: https://wiki.python.org/moin/HowTo/Sorting
def get_atm_id(item):
    return item.atm_id


def get_atm_key(item):
    return item.atm_key


def get_atm_grp(item):
    return item.atm_grp


def get_atm_sitnam(item):
    return item.sitnam


def get_atm_res(item):
    return item.res


def get_atm_chge(item):
    return item.chge


def get_atm_id1(item):
    return item.atm_id1


def get_atm_id2(item):
    return item.atm_id2


def get_atm_id3(item):
    return item.atm_id3


def get_atm_id4(item):
    return item.atm_id4
