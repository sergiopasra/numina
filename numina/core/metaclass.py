#
# Copyright 2008-2015 Universidad Complutense de Madrid
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

"""
Base metaclasses
"""

from .dataholders import Product
from .requirements import Requirement
import weakref


class StoreTypeAlt(type):
    """Metaclass for storing members."""
    def __new__(cls, classname, parents, attributes):

        newattr = {}
        __stored__ = weakref.WeakValueDictionary()
        for p in parents:
            stored = getattr(p, '__stored__', None)
            if stored:
                __stored__.update(stored)

        newattr.update(attributes)

        for name, val in newattr.items():
            if cls.exclude(name, val):
                nname, nval = cls.transform(name, val)
                __stored__[nname] = nval

        newattr['__stored__'] = __stored__

        return super(StoreTypeAlt, cls).__new__(cls, classname, parents, newattr)

    def __setattr__(self, name, value):
        # This is difficult to control safely...
        super(StoreTypeAlt, self).__setattr__(name, value)

    @classmethod
    def exclude(cls, name, value):
        return False

    @classmethod
    def transform(cls, name, value):
        return name, value


class StoreType(type):
    '''Metaclass for storing members.'''
    def __new__(cls, classname, parents, attributes):
        filter_out = {}
        filter_in = {}
        filter_in['__stored__'] = filter_out
        # Handle stored values from parents
        for p in parents:
            stored = getattr(p, '__stored__', None)
            if stored:
                filter_in['__stored__'].update(stored)

        for name, val in attributes.items():
            if cls.exclude(name, val):
                nname, nval = cls.transform(name, val)
                filter_out[nname] = nval
            else:
                filter_in[name] = val
        return super(StoreType, cls).__new__(
            cls, classname, parents, filter_in)

    def __setattr__(self, key, value):
        """Define __setattr__ in 'classes' created with this metaclass."""
        self._add_attr(key, value)

    def _add_attr(self, key, value):
        if self.exclude(key, value):
            self.__stored__[key] = value
        else:
            super(StoreType, self).__setattr__(key, value)

    @classmethod
    def exclude(cls, name, value):
        return False

    @classmethod
    def transform(cls, name, value):
        return name, value

class RecipeInOutType(StoreType):
    def __new__(cls, classname, parents, attributes):
        # Handle checkers defined in base class
        checkers = attributes.get('__checkers__', [])
        for p in parents:
            c = getattr(p, '__checkers__', [])
            checkers.extend(c)
        attributes['__checkers__'] = checkers
        obj = super(RecipeInOutType, cls).__new__(
            cls, classname, parents, attributes)
        return obj

    @classmethod
    def transform(cls, name, value):
        if value.dest is None:
            value.dest = name
        nname = value.dest
        return nname, value


class RecipeInputType(RecipeInOutType):
    """Metaclass for RecipeInput."""
    @classmethod
    def exclude(cls, name, value):
        return isinstance(value, Requirement)


class RecipeResultType(RecipeInOutType):
    """Metaclass for RecipeResult."""
    @classmethod
    def exclude(cls, name, value):
        return isinstance(value, Product)
