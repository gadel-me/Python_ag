#!/usr/bin/env python
from __future__ import print_function, division
import own_vmdfunctions as ovmd
import molecule

cbz = "/home/gadelmeier/Research/FORCE_FIELDS/AMBER/Molecules/CBZ_gaff2/4.4.VdW_vs_Coul/RESP/2H_esp_MP2_6-311Gsspp/RESP/CBZ-CBZ_2H.mol2"
molecule.load("mol2", cbz)
ovmd.label_charges(0, "charge")
