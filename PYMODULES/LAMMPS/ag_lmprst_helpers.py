#!/usr/bin/env python
from __future__ import print_function, division
import struct
import numpy as np


def import_lmprst(rst_lmp):
    f = open(rst_lmp)
    return f


def read_lmprst(rst_lmp, debug=False):
    """
    """
    #if debug:
    #    check_open(dcd_file)
    f = import_lmprst
    magic_string = f.read(len("LammpS RestartT"))
    endian = struct.unpack('i', f.read(4))
    version_numeric = struct.unpack('i', f.read(4))
    xperiodic = struct.unpack('i', f.read(4))

    # first 4 bytes of entry give length of whole block (bytes)
    first, = struct.unpack('i', f.read(4))
    return first

#write_string(VERSION,universe->version);
#write_int(SMALLINT,sizeof(smallint));
#write_int(IMAGEINT,sizeof(imageint));
#write_int(TAGINT,sizeof(tagint));
#write_int(BIGINT,sizeof(bigint));
#write_string(UNITS,update->unit_style);
#write_bigint(NTIMESTEP,update->ntimestep);
#write_int(DIMENSION,domain->dimension);
#write_int(NPROCS,nprocs);
#write_int_vec(PROCGRID,3,comm->procgrid);
#write_int(NEWTON_PAIR,force->newton_pair);
#write_int(NEWTON_BOND,force->newton_bond);
#write_int(XPERIODIC,domain->xperiodic);
#write_int(YPERIODIC,domain->yperiodic);
#write_int(ZPERIODIC,domain->zperiodic);
#write_int_vec(BOUNDARY,6,&domain->boundary[0][0]);
