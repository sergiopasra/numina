"""check_wlcalib_sp logs its information, instead of printing it with debugplot"""

import logging

import numpy as np

from numina.array.wavecalib.check_wlcalib import check_wlcalib_sp

LINES = np.array([5050.0, 5130.0, 5210.5, 5333.0, 5401.0, 5520.0, 5640.0, 5777.0, 5850.0, 5930.0])


def arc_spectrum():
    wl = 5000.0 + np.arange(1000) * 1.0
    return 10 + sum(1000 * np.exp(-0.5 * ((wl - line) / 1.5) ** 2) for line in LINES)


def test_messages_are_logged(caplog, capsys):
    caplog.set_level(logging.DEBUG, logger="numina.array.wavecalib.check_wlcalib")
    check_wlcalib_sp(arc_spectrum(), crpix1=1.0, crval1=5000.2, cdelt1=1.0, wv_master=LINES, debugplot=0)
    assert ">>> Number of identified lines........: 10" in caplog.text
    assert capsys.readouterr().out == ""
