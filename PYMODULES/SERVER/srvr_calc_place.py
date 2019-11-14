#!/usr/bin/env python

import os
import socket

"""
Search on which server we are currently.
"""

__version__ = "2017-05-30"


def get_fasttmp():
    """
    Sort out if a $FASTTMP is existing and its location
    """
    hostname = socket.gethostname()

    try:
        fasttmp = os.environ["FASTTMP"]  # lima or emmy
    except KeyError:
        if hostname in ("wood3" or "zahnserver"):
            fasttmp = "/scratch/"  # fasttmp = scratch folder
        else:
            fasttmp = None  # local host

    return fasttmp
