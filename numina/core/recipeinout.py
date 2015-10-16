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
Recipe Input
"""

from six import with_metaclass
import yaml

from .metaclass import RecipeInputType, RecipeResultType
import numina.store.dump


class RecipeInOut(object):

    def __init__(self, *args, **kwds):
        super(RecipeInOut, self).__init__()
        # Used to hold set values
        self._numina_desc_val = {}

        for key, val in kwds.items():
            setattr(self, key, kwds[key])

        self._finalize()

    def _finalize(self):
        """Access all the instance descriptors

        This wil trigger an exception if a required
        parameter is not set
        """
        all_msg_errors = []

        for key in self.stored():
            try:
                getattr(self, key)
            except ValueError as err:
                all_msg_errors.append(err.message)

        # Raises a list of all the missing entries
        if all_msg_errors:
            raise ValueError(all_msg_errors)

    def attrs(self):
        return self._numina_desc_val

    @classmethod
    def stored(cls):
        return cls.__numina_stored__

    def validate(self):
        """Validate myself."""

        for key, req in self.stored().items():
            val = getattr(self, key)
            req.validate(val)

        # Run checks defined in __checkers__
        self._run_checks()

    def _run_checks(self):
        checkers = getattr(self, '__checkers__', [])

        for check in checkers:
            check.check(self)


class RecipeInput(with_metaclass(RecipeInputType, RecipeInOut)):
    """RecipeInput base class"""
    pass


class BaseRecipeResult(object):
    def __new__(cls, *args, **kwds):
        return super(BaseRecipeResult, cls).__new__(cls)

    def __init__(self, *args, **kwds):
        super(BaseRecipeResult, self).__init__()


class ErrorRecipeResult(BaseRecipeResult):
    def __init__(self, errortype, message, traceback):
        super(ErrorRecipeResult, self).__init__()
        self.errortype = errortype
        self.message = message
        self.traceback = traceback

    def __repr__(self):
        sclass = type(self).__name__
        fmt = "%s(errortype=%r, message='%s')"
        return fmt % (sclass, self.errortype, self.message)

    def store_to(self, where):
        with open(where.result, 'w+') as fd:
            yaml.dump(where.result, fd)

        return where.result


class RecipeResult(with_metaclass(RecipeResultType, RecipeInOut, BaseRecipeResult)):

    def __repr__(self):
        sclass = type(self).__name__
        full = []
        for key, val in self.stored().items():
            full.append('%s=%r' % (key, val))
        return '%s(%s)' % (sclass, ', '.join(full))

    def store_to(self, where):

        saveres = {}
        for key, prod in self.stored().items():
            val = getattr(self, key)
            where.destination = prod.dest
            saveres[key] = numina.store.dump(prod.type, val, where)

        with open(where.result, 'w+') as fd:
            yaml.dump(saveres, fd)

        return where.result


class define_result(object):
    """Recipe decorator."""
    def __init__(self, resultClass):
        if not issubclass(resultClass, RecipeResult):
            msg = '%r does not derive from RecipeResult' % resultClass
            raise TypeError(msg)
        self.klass = resultClass

    def __call__(self, klass):
        klass.RecipeResult = self.klass
        return klass


class define_input(object):
    """Recipe decorator."""
    def __init__(self, inputClass):
        if not issubclass(inputClass, RecipeInput):
            fmt = '%r does not derive from RecipeInput'
            msg = fmt % inputClass
            raise TypeError(msg)
        self.klass = inputClass

    def __call__(self, klass):
        klass.RecipeInput = self.klass
        return klass


define_requirements = define_input


class add_requirement(object):
    def __init__(self, **kwds):

        self.ext = {}
        for key, val in kwds.items():
            # FIXME validate these inputs
            self.ext[key] = val

    def __call__(self, klass):
        Class = klass.RecipeInput
        for key, val in self.ext.items():
            setattr(Class, key, val)
        return klass


class add_product(object):
    def __init__(self, **kwds):

        self.ext = {}
        for key, val in kwds.items():
            # FIXME validate these inputs
            self.ext[key] = val

    def __call__(self, klass):
        Class = klass.RecipeResult
        for key, val in self.ext.items():
            setattr(Class, key, val)
        return klass