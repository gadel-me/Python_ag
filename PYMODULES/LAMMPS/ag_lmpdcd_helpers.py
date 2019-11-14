#!/usr/bin/env python

import struct
import numpy as np


def read_record(dcd_file, debug=False):
    """
    First and list integer of each record contains its length in byte.
    (First and last integer (8 bytes) are excluded!)
    e.g.  84     CORD 1 2 4 5 3 5      84
         187    2 395 9 59 82 87 2    187
           4      298985                4
    """
    if debug:
        check_open(dcd_file)

    # first 4 bytes of entry give length of whole block (bytes)
    first, = struct.unpack('i', dcd_file.read(4))
    # actual record (of 'first' bytes length)
    record = dcd_file.read(first)
    # last 4 bytes signal end of entry
    last, = struct.unpack('i', dcd_file.read(4))

    if debug:
        print("Length of block in bytes: ", first)

    # values of first and last entry should be the same
    if first != last:
        raise IOError("First ({}) and last ({}) 4 bytes are not the same!".format(
            first, last))

    return record


def check_open(dcd_file):
    if dcd_file.closed:
        raise ValueError("I/O operation on closed file")


def deploy_array(num_frames, natoms):
    """
    Preallocate memory for upcoming arrays
    """
    x_coordset = np.zeros((num_frames, natoms))
    y_coordset = np.zeros((num_frames, natoms))
    z_coordset = np.zeros((num_frames, natoms))

    return (x_coordset, y_coordset, z_coordset)


def reshape_arguments(sframe, nframes, step, frame_start, frame_stop, key):
    """
    Helper function for method 'read_frames' of class DCDFile (ag_lmpdcd-module).

    Convert frame_start and frame_stop so that they may contain values such as
    -1 or None, the list indices or the actual frame numbers.
    """
    # /// CHECK INPUT
    if frame_start == frame_stop:
        raise ValueError("Arguments 'frame' and 'to_frame' must not be " +
                         "the same value (last is always exclusive)!")

    # convert frame numbers to indices (subtract number of first DCD frame)
    if key == "step":
        nframes *= step  # number of total steps
    else:
        sframe = 0  # first frame has index 0
        step = 1  # each index has step size of 1

    # make values < 0 and None available for frame_start
    if frame_start is None:
        fsta = sframe
    elif frame_start < 0:
        # start counting from behind (-1 = last frame, -2 = second last)
        fsta = (nframes+step) + (frame_start*step)
    else:
        fsta = frame_start

    # make values < 0 and None available for frame_start
    if frame_stop is None:
        fsto = nframes + step  # add one step to be exclusive for last frame in read_frames-method
    elif frame_stop < 0:
        # make it possible to count from behind (-1 = last frame, -2 = second last)
        fsto = (nframes+step) + (frame_stop*step)
    elif frame_stop == "no":
        fsto = fsta + step
    else:
        fsto = frame_stop

    # check if indices are o.k.
    check_for_err(sframe, nframes, step, fsta)
    if frame_stop is None:
        nframes += step  # add one frame (read_frames-method reasons)
    check_for_err(sframe, nframes, step, fsto)

    # convert values to indices
    if frame_start < 0:
        fsta = int(fsta/step)
    else:
        fsta = int((frame_start-sframe)/step)

    if frame_stop < 0:
        fsto = int(fsto/step)
    else:
        fsto  = int((fsto-sframe)/step)

    return (fsta, fsto)


def check_for_err(sf, nfr, ste, fst):
    """
    sf:     sframe
    nfr:    nframes
    ste:    step
    fst:   fstart/fstop
    """

    if not (sf <= fst <= nfr):
        raise ValueError("{1} not between {0} and {2}! ".format(
                         sf, fst, nfr) +
                         "Choose another frame as starting point.")
