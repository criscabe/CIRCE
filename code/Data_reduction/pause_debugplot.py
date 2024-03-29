#
# Copyright 2019-2022 Universidad Complutense de Madrid
#
# This file software has been employed to reduce infrared 
# raw data from the CIRCE camera at GTC (see Cabello et al. 2022)
#
# Authors: Cristina Cabello (criscabe@ucm.es), 
#          Nicolás Cardiel (cardiel@ucm.es)
#          Jesús Gallego (j.gallego@ucm.es)
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#
from __future__ import division
from __future__ import print_function


def pause_debugplot(debugplot):
    """Ask the user to press RETURN to continue after plotting.

    Parameters
    ----------
    debugplot : int
        Determines whether intermediate computations and/or plots
        are displayed:
        00 : no debug, no plots
        01 : no debug, plots without pauses
        02 : no debug, plots with pauses
        10 : debug, no plots
        11 : debug, plots without pauses
        12 : debug, plots with pauses

    """

    if debugplot == 2 or debugplot == 12:
        try:
            input("\nPress RETURN to continue...")
        except SyntaxError:
            pass

        print(' ')
     
        
