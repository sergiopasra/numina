#
# Copyright 2008-2017 Universidad Complutense de Madrid
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

import inspect

from numina.util.parser import parse_arg_line
from numina.exceptions import NoResultFound
from numina.datamodel import DataModel


class DataTypeBase(object):
    """Base class for input/output types of recipes.

    """
    def __init__(self, *args, **kwds):
        if 'datamodel' in kwds:
            datamodel = kwds['datamodel']
            if inspect.isclass(datamodel):
                self.datamodel = datamodel()
            else:
                self.datamodel = datamodel
        else:
            self.datamodel = DataModel()

    def __getstate__(self):
        return {}

    def __setstate__(self, state):
        pass

    def query(self, name, dal, ob, options=None):

        try:
            return self.query_on_ob(name, ob)
        except NoResultFound:
            pass

        #param = dal.search_param_req_tags(req, ob.instrument,
        #                                      ob.mode, ob.tags, ob.pipeline)
        param = dal.search_parameter(name, self, ob)
        return param.content


    def query_on_ob(self, key, ob):
        # First check if the requirement is embedded
        # in the observation result
        # it can happen in GTC

        try:
            return getattr(ob, key)
        except AttributeError:
            raise NoResultFound("DataType.query_on_ob")

    def on_query_not_found(self, notfound):
        pass

    @classmethod
    def isproduct(cls):
        """Check if the DataType is the product of a Recipe"""
        return False

    def __repr__(self):
        sclass = type(self).__name__
        return "%s()" % (sclass, )

    def name(self):
        """Unique name of the datatype"""
        return self.__repr__()

    @classmethod
    def from_name(cls, name):
        # name is in the form Class(arg1=val1, arg2=val2)
        # find first (
        fp = name.find('(')
        sp = name.find(')')
        if (fp == -1) or (sp == -1) or (sp < fp):
            # paren not found
            klass = name
            kwargs = {}
        else:
            # parse things between parens
            klass = name[:fp]
            fargs = name[fp+1:sp]
            kwargs = parse_arg_line(fargs)

        if klass == cls.__name__:
            # create thing
            obj = cls.__new__(cls)
            obj.__init__(**kwargs)

            return obj
        else:
            raise TypeError(name)

    def extract_tags(self, obj):
        """Extract tags from serialized file"""
        meta_info = self.extract_meta_info(obj)
        return meta_info['tags']

    @staticmethod
    def create_meta_info():
        """Create metadata structure"""
        result = {}
        result['instrument'] = ''
        result['uuid'] = ''
        result['tags'] = {}
        result['type'] = ''
        result['mode'] = ''
        result['observation_date'] = ""
        result['origin'] = {}
        return result

    def extract_meta_info(self, obj):
        """Extract metadata from serialized file"""
        result = self.create_meta_info()
        return result