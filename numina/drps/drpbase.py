#
# Copyright 2011-2016 Universidad Complutense de Madrid
#
# This file is part of Numina
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""DRP storage class"""

import warnings


class DrpBase(object):
    """Store DRPs, base class"""

    def __init__(self):
        pass

    def query_by_name(self, name):
        return None

    def query_all(self):
        return {}

    @staticmethod
    def instrumentdrp_check(drpins, entryname):
        import numina.core.pipeline
        if isinstance(drpins, numina.core.pipeline.InstrumentDRP):
            if drpins.name == entryname:
                return True
            else:
                msg = 'Entry name "{}" and DRP name "{}" differ'.format(entryname, drpins.name)
                warnings.warn(msg, RuntimeWarning)
                return False
        else:
            msg = 'Object {0!r} does not contain a valid DRP'.format(drpins)
            warnings.warn(msg, RuntimeWarning)
            return False


class DrpGeneric(DrpBase):
    """Store DRPs.in a dictionary"""

    def __init__(self, drps=None):
        super(DrpGeneric, self).__init__()
        self.drps = drps

    def query_by_name(self, name):
        """Query DRPs in internal storage by name

        """
        return self._drps[name]

    def query_all(self):
        """Return all available DRPs in internal storage"""

        return self._drps

    @property
    def drps(self):
        return self._drps

    @drps.setter
    def drps(self, newdrps):
        self._drps = {} if newdrps is None else newdrps