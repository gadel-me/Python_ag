from __future__ import print_function, division
import numpy as np
#import copy
import md_stars as mds
import md_box as mdb
import md_universe as mdu
import md_elements as mde

__version__ = "2018-04-17"


class XYZ(mdu.Universe):
    """
    Read/Write XYZ-Files.
    Sources:    https://en.wikipedia.org/wiki/XYZ_file_format
                http://openbabel.org/wiki/XYZ_%28format%29
    """
    def __init__(self):
        """
        Load universe first.
        """
        mdu.Universe.__init__(self)

    def read_xyz(self, xyzfile, overwrite_data=False):
        """
        Read the contents of a xyz-file.
        XYZ-Format: https://en.wikipedia.org/wiki/XYZ_file_format
        """
        with open(xyzfile, "r") as xyz_in:
            for line in xyz_in:
                num_atms = int(line.split()[0])  # line with number of atoms (mandatory)
                comment_line = xyz_in.next()  # comment line (may be empty)

                # read box information if there is any
                if "Boxtype" in comment_line:
                    comment_line = comment_line.split()
                    box = [None if i == "None" else None if round(float(i), 8) == 0 else float(i) for i in comment_line[3::2]]
                    print(box)

                    if comment_line[1] == "cartesian":
                        current_box = mdb.Box(boxtype=comment_line[1],
                                              crt_a=box[0],
                                              crt_b=box[1],
                                              crt_c=box[2])
                    elif comment_line[1] == "lattice":
                        current_box = mdb.Box(boxtype=comment_line[1],
                                              ltc_a=box[0],
                                              ltc_b=box[1],
                                              ltc_c=box[2],
                                              ltc_alpha=box[3],
                                              ltc_beta=box[4],
                                              ltc_gamma=box[5])
                    else:
                        current_box = mdb.Box(boxtype=comment_line[1],
                                              lmp_xlo=box[0],
                                              lmp_xhi=box[1],
                                              lmp_ylo=box[2],
                                              lmp_yhi=box[3],
                                              lmp_zlo=box[4],
                                              lmp_zhi=box[5],
                                              lmp_xy=box[6],
                                              lmp_xz=box[7],
                                              lmp_yz=box[8])

                    del (box, comment_line)

                    # overwrite given boxes
                    if overwrite_data is True:
                        self.ts_boxes = []

                    self.ts_boxes.append(current_box)

                cframe = []

                # parse coordinates section
                for iid in xrange(num_atms):
                    catm = xyz_in.next().split()
                    csitnam = catm[0]
                    ccoords = np.array([float(i) for i in catm[1:5]])

                    # charge of current atom
                    ccharge = None
                    if len(catm) > 4:
                        ccharge = float(catm[4])

                    # check if an instance of Atom with index iid already exists;
                    # overwrite info if it does or create a new one if it does not
                    try:
                        self.atoms[iid]

                        # overwrite data
                        if overwrite_data is True:
                            self.atoms[iid].atm_id = iid
                            self.atoms[iid].sitnam = csitnam

                            if ccharge is not None:
                                self.atoms[iid].chge = ccharge

                        # complement data
                        else:

                            if not hasattr(self.atoms[iid], "atm_id"):
                                self.atoms[iid].atm_id = iid

                            if not hasattr(self.atoms[iid], "sitnam"):
                                self.atoms[iid].sitnam = csitnam

                            if not hasattr(self.atoms[iid], "chge") and ccharge is not None:
                                self.atoms[iid].chge = ccharge

                    except IndexError:
                            catm = mds.Atom(atm_id=iid,
                                            sitnam=csitnam)
                            if ccharge is not None:
                                catm.chge = ccharge
                            self.atoms.append(catm)

                    cframe.append(ccoords)

                # append current frame
                self.ts_coords.append(cframe)

    def write_xyz(self, xyz_file, title, guess_element=False, *frame_ids):
        """
        Write a xyz-file.
        """
        if frame_ids == ():
            frame_ids = [0]

        num_atms = len(self.atoms)

        with open(xyz_file, "w") as xyz_out:

            for frame_id in frame_ids:
                xyz_out.write("{}\n".format(num_atms))
                xyz_out.write("Step: {} - {}".format(frame_id, title))
                xyz_out.write("\n")

                # write box info to comment line
                if self.ts_boxes != [] and self.ts_boxes is not None:
                    xyz_out.write("Boxtype {} ".format(self.ts_boxes[frame_id].boxtype))

                    if self.ts_boxes[frame_id].boxtype == "cartesian":
                        xyz_out.write("a: {} b: {} c: {}\n".format(self.ts_boxes[frame_id].crt_a,
                                                                   self.ts_boxes[frame_id].crt_b,
                                                                   self.ts_boxes[frame_id].crt_c))
                    elif self.ts_boxes[frame_id].boxtype == "lattice":
                        xyz_out.write(("a: {} b: {} c: {} "
                                       "alpha: {} beta: {} gamma: {}\n").format(self.ts_boxes[frame_id].ltc_a,
                                                                                self.ts_boxes[frame_id].ltc_b,
                                                                                self.ts_boxes[frame_id].ltc_c,
                                                                                self.ts_boxes[frame_id].ltc_alpha,
                                                                                self.ts_boxes[frame_id].ltc_beta,
                                                                                self.ts_boxes[frame_id].ltc_gamma))
                    else:
                        xyz_out.write(("xlo {} xhi {} "
                                       "ylo {} yhi {} "
                                       "zlo {} zhi {} "
                                       "xy {} xz {} yz {}\n").format(self.ts_boxes[frame_id].lmp_xlo,
                                                                     self.ts_boxes[frame_id].lmp_xhi,
                                                                     self.ts_boxes[frame_id].lmp_ylo,
                                                                     self.ts_boxes[frame_id].lmp_yhi,
                                                                     self.ts_boxes[frame_id].lmp_zlo,
                                                                     self.ts_boxes[frame_id].lmp_zhi,
                                                                     self.ts_boxes[frame_id].lmp_xy,
                                                                     self.ts_boxes[frame_id].lmp_xz,
                                                                     self.ts_boxes[frame_id].lmp_yz))

                for catm, ccoords in zip(self.atoms, self.ts_coords[frame_id]):
                    # write atom name
                    if hasattr(catm, "sitnam") is True and guess_element is False:
                        xyz_out.write("{} ".format(catm.sitnam))
                    else:
                        # guess element by mass
                        cmass = round(self.atm_types[catm.atm_key].weigh, 1)

                        try:
                            sitnam = mde.elements[cmass]
                        except AttributeError:
                            sitnam = "X"

                        xyz_out.write("{} ".format(sitnam))

                    # write coordinates
                    xyz_out.write("{:> 10.5f} {:> 10.5f} {:> 10.5f} ".format(
                        ccoords[0], ccoords[1], ccoords[2]))

                    # write charges as additional info if they were given
                    if hasattr(catm, "chge") is True:
                        xyz_out.write("{:> 10.5f} ".format(catm.chge))

                    # write a newline after current entry is done
                    xyz_out.write("\n")
