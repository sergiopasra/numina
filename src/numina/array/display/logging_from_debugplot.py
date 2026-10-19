#
# Copyright 2018-2023 Universidad Complutense de Madrid
#
# This file is part of Numina
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

import logging


def logging_from_debugplot(debugplot, fmt=None):
    """Set debugging level based on debugplot value.

    Parameters
    ----------
    debugplot : int
            Debugging level for messages and plots. For details see
            'numina.array.display.pause_debugplot.py'.
    fmt : str, optional
            Format of the messages, the default of logging.basicConfig
            if None. With '%(message)s' the messages are shown as with print.
    """

    kwargs = {} if fmt is None else {"format": fmt}
    if isinstance(debugplot, int):
        if abs(debugplot) >= 10:
            logging.basicConfig(level=logging.DEBUG, **kwargs)
        else:
            logging.basicConfig(level=logging.INFO, **kwargs)
        # the messages of matplotlib, as finding fonts, are not of interest
        logging.getLogger("matplotlib").setLevel(logging.ERROR)
    else:
        raise ValueError("Unexpected debugplot=" + str(debugplot))
