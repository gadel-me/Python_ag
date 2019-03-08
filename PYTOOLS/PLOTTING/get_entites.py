from __future__ import print_function, division
import pdb
import argparse
import numpy as np
import ag_geometry as agm
import ag_clmsv
import ag_unify_md as agum


def chunks(l, n):
    """
    Yield successive n-sized chunks from l.

    Sources : https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
    """
    return [l[i:i + n] for i in xrange(0, len(l), n)]


def get_entity(crds, idxs, molecule_offset):
    """
    Get measured entity from given coordinates.

    Parameters
    ----------
    crds : numpy.ndarray
        coordinates to measure entity from
    idxs : list
        list of atom indices which form the entity
    molecule_offset : int
        number of atoms of a molecule/other structure which increments the atom
        indices

    Returns
    -------
    measure_entites : list
        all measured entities from the given coordinates

    """
    measured_entites = []
    len_coords = len(crds)

    while idxs[-1] < len_coords:

        # Bond
        if len(idxs) == 2:
            measured_entity = np.linalg.norm(
                crds[idxs[0]] - crds[idxs[1]])
        # angle
        elif len(idxs) == 3:
            measured_entity = np.degrees(agm.get_angle(
                crds[idxs[0]], crds[idxs[1]],
                crds[idxs[2]]))
        # dihedral or improper
        elif len(idxs) == 4:
            measured_entity = np.degrees(agm.get_dihedral(
                crds[idxs[0]], crds[idxs[1]],
                crds[idxs[2]], crds[idxs[3]]))
        else:
            raise Warning("***Too many atom indices! Unclear what to measure")

        measured_entites.append(measured_entity)
        # increment idxs by given offset
        idxs = [i + molecule_offset for i in idxs]

    return measured_entites


def process_entity_idxs(data, crds, entity_idxs, entity_type, increment=0):
    """
    Get all entities from a given list of entity indices.

    Parameters
    ----------
    data : list
        list to append the data of each entity to
    crds : np.ndarray
        coordinates to process

    entity_idxs : list
        list of lists with entity_idxs, e.g. [[1,2], [3,4], [4,5],...]

    entity_type : str
        type of entity to create the name from;
        only 'bonds', 'angles' and 'dihedrals' are allowed

    """
    complete_entity = {}

    if entity_type.lower() == "bonds":
        entity_name = "{c[0]} {c[1]}"
    elif entity_type.lower() == "angles":
        entity_name = "{c[0]} {c[1]} {c[2]}"
    else:
        entity_name = "{c[0]} {c[1]} {c[2]} {c[3]}"

    for idxs in entity_idxs:
        # get all values for current entity with the currently given indices
        entity_vals = get_entity(crds, idxs, ARGS.offset)

        # current entity to dict
        entry_name = entity_name.format(c=[i + increment for i in idxs])
        complete_entity[entry_name] = entity_vals

    data.append(complete_entity)


if __name__ == "__main__":
    #=== parse arguments ===
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("-lmpdat", default=None,
                        help="Lammps' data-file")
    PARSER.add_argument("-dcd",
                        help="dcd file from minimization run of the cbz dimer system")
    PARSER.add_argument("-bonds", nargs="*", type=int, default=None,
                        help="indices of atoms that form bonds")
    PARSER.add_argument("-angles", nargs="*", type=int, default=None,
                        help="indices of atoms that form angles")
    PARSER.add_argument("-dihedrals", nargs="*", type=int, default=None,
                        help="indices of atoms that form dihedrals/impropers")
    PARSER.add_argument("-offset", type=int, default=30,
                        help="number of atoms a molecule/scan-block consists of")

    PARSER.add_argument("-o", default="DEFAULTNAME", help="Prefix of output files")
    ARGS = PARSER.parse_args()

    # container
    COMPLETE_DATA = ag_clmsv.Clmsv()
    MOLSYS = agum.Unification()

    if ARGS.lmpdat is not None:
        MOLSYS.read_lmpdat(ARGS.lmpdat)

    if ARGS.dcd is not None:
        # read last frame from dcd file
        MOLSYS.import_dcd(ARGS.dcd)
        MOLSYS.read_frames(frame=-2)

    #=== split the indices ===
    if ARGS.bonds is not None:
        ENTITY_IDXS = chunks(ARGS.bonds, 2)
        process_entity_idxs(COMPLETE_DATA.data,
                            MOLSYS.ts_coords[-1],
                            ENTITY_IDXS, "bonds", increment=1)

    if ARGS.angles is not None:
        ENTITY_IDXS = chunks(ARGS.angles, 3)
        process_entity_idxs(COMPLETE_DATA.data,
                            MOLSYS.ts_coords[-1],
                            ENTITY_IDXS, "angles", increment=1)

    if ARGS.dihedrals is not None:
        ENTITY_IDXS = chunks(ARGS.dihedrals, 4)
        process_entity_idxs(COMPLETE_DATA.data,
                            MOLSYS.ts_coords[-1],
                            ENTITY_IDXS, "dihedrals", increment=1)

    COMPLETE_DATA.write_clmsv(ARGS.o)
