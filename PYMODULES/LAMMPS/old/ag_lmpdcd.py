#!/usr/bin/env python
from __future__ import print_function, division
import struct
import numpy as np
import md_box as mdb
import md_universe as mdu
import ag_lmpdcd_helpers as agldh

__version__ = "2017-04-27"


class LmpDCD(mdu.Universe):
    """
    Read and write the "DCD" binary trajectory file format used by LAMMPS.
    Other formats are not considered (yet).
    Sources:    http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
                https://docs.python.org/2/library/struct.html#functions-and-exceptions
                http://www.devdungeon.com/content/working-binary-data-python
                http://stackoverflow.com/questions/1035340/reading-binary-file-in-python-and-looping-over-each-byte/1035360#1035360
                http://stackoverflow.com/questions/38297929/why-does-unpacking-a-struct-result-in-a-tuple
                http://prody.csb.pitt.edu/_modules/prody/trajectory/dcdfile.html#codemodal
    #a = LmpDCD(dcd)
    #print("1,2")
    #a.read_frames(1, 2)
    #print("1,None")
    #a.read_frames(1, None)
    #print("2, None")
    #a.read_frames(2, None)
    #print("None, 10")
    #a.read_frames(None, 10)
    #print("None, -1")
    #a.read_frames(None, -1)
    #print("2001 2001")
    #a.read_frames(frame=0, to_frame=None, frame_by="index")
    #a.jump_to_first_frame()
    #a.read_frames(frame=-505, to_frame=-501, frame_by="index")
    #a.jump_to_first_frame()
    #a.read_frames(frame=820000, to_frame=825000, frame_by="step", verbose=True)
    #a.jump_to_first_frame()
    #a.read_frames(frame=0, to_frame="no", frame_by="index", verbose=True)
    #a.jump_to_first_frame()
    #a.read_frames(frame=0, to_frame=300, frame_by="index", verbose=True)
    ##a.jump_to_first_frame()
    #a.read_frames(frame=400, to_frame=500, frame_by="index", verbose=True)
    #a.jump_to_first_frame()
    #a.append_dcds(dcd2, dcd, frame=0, to_frame=5, frame_by="index", verbose=True)
    #a.read_frames(0, None)

    #b = LmpDCD()
    #b.import_dcd(dcd)
    #b.read_frames()

    """
    def __init__(self):
        """
        Lorem ipsum...
        """
        mdu.Universe.__init__(self)

    def import_dcd(self, dcd):
        """
        Open DCD and read the remarks. This function must be called before
        anything else that has to do anything with file reading.
        """
        self._dcdfile = open(dcd, 'rb')
        self._read_header()
        self._read_title()
        # get position in file (bytes?) after header and title
        self._pos_1 = self._dcdfile.tell()

    def jump_to_first_frame(self):
        """
        Rewind to the first frame
        """
        self._dcdfile.seek(self._pos_1)

    def _read_header(self):
        """
        Read header block.
        """
        self.extra_blck = None  # only if unit cell
        self.has_4dims  = None  # purpose unknown
        self.is_charmm  = False

        # header block
        hdr_blck = agldh.read_record(self._dcdfile)
        hdr = struct.unpack('4c9if10i', hdr_blck)
        self.nframes = hdr[4]  # total number of frames
        self.sframe  = hdr[5]  # number of start frame
        self.step    = hdr[6]  # number of frames between each frame
        self.lframe  = hdr[7]  # number of last frame

        # check charmm-formatting
        if hdr[23] != 0:
            self.is_charmm = True
            self.extra_blck = hdr[14]
            # some dcd files have 4 dimensions?
            self.has_4dims = hdr[15]
            print("***Reading Charmm formatted DCD with 4 dimensions!***")

    def _read_title(self):
        """
        Read title blocks (only possible if header block was read before!)
        """
        # 2 title lines
        title_1 = agldh.read_record(self._dcdfile)  # 1st title block
        title_2 = agldh.read_record(self._dcdfile)  # 2nd title block
        self.natoms, = struct.unpack("i", title_2)
        print("   Remark 1: {}\n   Remark 2: Number of Atoms: {}".format(
              title_1, self.natoms))

    def _read_frame(self):
        """
        Read one frame.
        Layout of unitcell is [A, alpha, B, beta, gamma, C] (Historical reasons)
        """
        # read cell information
        if self.extra_blck:
            cur_blck = agldh.read_record(self._dcdfile)
            cur_cell = struct.unpack("6d", cur_blck)

        # read coordinate-sets
        x_coordset = agldh.read_record(self._dcdfile)
        y_coordset = agldh.read_record(self._dcdfile)
        z_coordset = agldh.read_record(self._dcdfile)
        x = np.fromstring(x_coordset, dtype=np.dtype('f'), count=self.natoms)
        y = np.fromstring(y_coordset, dtype=np.dtype('f'), count=self.natoms)
        z = np.fromstring(z_coordset, dtype=np.dtype('f'), count=self.natoms)

        # 4th dimension given? (has also to be read)
        if self.has_4dims:
            agldh.read_record(self._dcdfile)
            #dims_4_blck = agldh.read_record(self._dcdfile)
            #dunno = struct.unpack(str(self.natoms)+"i", dims_4_blck)

        return(x, y, z, cur_cell)

    def _skip_frame(self):
        if self.extra_blck:
            agldh.read_record(self._dcdfile)

        agldh.read_record(self._dcdfile)
        agldh.read_record(self._dcdfile)
        agldh.read_record(self._dcdfile)

        if self.has_4dims:
            agldh.read_record(self._dcdfile)

    def read_frames(self, frame=None, to_frame=-1, frame_by="index", verbose=False):
        """
        Sources:    https://github.com/MDAnalysis/mdanalysis/issues/187
        """
        if verbose:
            print("***Verbose: MB already read: {:.2f} MiB".format(self._dcdfile.tell()/1000000))

        M_PI_2 = np.pi/2
        # convert input to corresponding indices
        frm, to_frm = agldh.reshape_arguments(self.sframe, self.nframes,
                                              self.step, frame, to_frame,
                                              frame_by)

        # preallocate memory for arrays to come
        if to_frame is None:
            num_frames = abs(frm - to_frm) - 1
        else:
            num_frames = abs(frm - to_frm)  # one frame is always read

        print("***Info: Reading: Frame (start): {}, ToFrame (excluded): {}, NumFrames: {}".format(frm, to_frm, num_frames))
        x_set, y_set, z_set = agldh.deploy_array(num_frames, self.natoms)
        # fill arrays with coordinates
        ptr = 0  # pointer to place data in right position of array

        for frame_num in xrange(self.nframes):

            if frm <= frame_num < to_frm:

                # coords
                x, y, z, cur_box = self._read_frame()
                x_set[ptr] = x
                y_set[ptr] = y
                z_set[ptr] = z

                # create box and append to other boxes of trajectory
                #TODO: Check if angles are right this way with triclinic cell
                alpha = np.radians(90.0 - np.arcsin(cur_box[4])*90.0/M_PI_2)  # cosAB
                beta  = np.radians(90.0 - np.arcsin(cur_box[3])*90.0/M_PI_2)  # cosAC
                gamma = np.radians(90.0 - np.arcsin(cur_box[1])*90.0/M_PI_2)  # cosBC
                cur_box = mdb.Box(ltc_alpha=alpha,
                                  ltc_beta=beta,
                                  ltc_gamma=gamma,
                                  ltc_a=cur_box[0],
                                  ltc_b=cur_box[2],
                                  ltc_c=cur_box[5],
                                  boxtype="lattice")

                # convert lattice-box-type to lammps-box-type
                cur_box.box_lat2lmp()
                self.ts_boxes.append(cur_box)

                # get step numbers per frame so we can later access them if wanted
                #self.ts_steps.append(frame_num*self.step+self.sframe)
                ptr += 1
            elif frame_num > to_frm:
                break
            else:
                self._skip_frame()  # keep reading until first is reached

        # concatenate all arrays to (1,3)-arrays
        coordinates = np.stack((x_set, y_set, z_set), axis=-1)

        # append coordinates to universe ts-coordinates
        for i in coordinates:
            self.ts_coords.append(i)

    def append_dcds(self, *dcd_files, **read_frame_args):
        """
        Read and append another DCD to an existing one/s.
        Only dcd may appended, if the number of atoms is the same!

        dcds:   str; dcd file name(s) to append
        read_frame_args: same as for method read_frames
        """
        cur_num_atoms = self.natoms  # same atoms as in data

        for cur_dcd in dcd_files:
            c_dcd = LmpDCD()
            c_dcd.import_dcd(cur_dcd)

            # check if number of atoms is the same as in data/previous dcd file
            if cur_num_atoms != c_dcd.natoms:
                raise RuntimeError("Different number of atoms in DCD-files!")

            cur_num_atoms = c_dcd.natoms
            c_dcd.read_frames(**read_frame_args)  # use same kwargs as append_dcds
            do_not_append = []

            # find duplicate step-entries, save indices
            for k, i in enumerate(c_dcd.ts_steps):
                for j in self.ts_steps:
                    if i == j:
                        do_not_append.append(k)

            # append ts-steps, ts-boxes and ts-coordinates to universe
            for iidx, istp in enumerate(c_dcd.ts_steps):
                if iidx not in do_not_append:
                    self.ts_steps.append(istp)
                    self.ts_boxes.append(c_dcd.ts_boxes[iidx])
                    self.ts_coords.append(c_dcd.ts_coords[iidx])

            c_dcd.close_dcd()

    def close_dcd(self):
        """
        Close dcdfile if still open.
        """
        if not self._dcdfile.closed:
            print("***Info: Closing file: {}.".format(self._dcdfile))
            self._dcdfile.close()

    def write_dcd(self, *frames):
        """
        Write frames to DCD-file. TBD.
        """
        pass
