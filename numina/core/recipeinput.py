#
# Copyright 2008-2014 Universidad Complutense de Madrid
#
# This file is part of Numina
#
# Numina is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Numina is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Numina.  If not, see <http://www.gnu.org/licenses/>.
#

from numina.core import RequirementParser
from numina.core.load import dict_requirement_lookup


class RecipeInputBuilder(object):

    def build(self, workenv, klass, mreqs):
        lc = dict_requirement_lookup(mreqs)
        rp = RequirementParser(klass, lookup=lc)
        requires = rp.parse()

        workenv.sane_work()
        workenv.copyfiles(mreqs['obresult'], requires)

        return requires