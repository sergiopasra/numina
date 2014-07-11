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


import inspect

from .metaclass import MapStoreType
from .products import DataProduct, QualityControlProduct
from .types import NullType, PlainPythonType
from .types import ListOf

class Product(object):
    '''Product holder for RecipeResult.'''
    def __init__(self, product_type, description='', validate=False, 
            dest=None, optional=False, default=None, *args, **kwds):

        self.type = product_type
        if isinstance(self.type, Optional):
            self.type = self.type.product_type
            self.optional = True

        if self.type is None:
            self.type = NullType()
        elif self.type in [bool, str, int, float, complex]:
            self.type = PlainPythonType(ref=product_type())
        elif isinstance(self.type, ListOf):
            self.type = self.type
        else:
            if inspect.isclass(self.type):
                self.type = self.type()
                
            if isinstance(self.type, DataProduct):
                pass
            else:
                raise TypeError('product_type must be of class DataProduct')
        
        self.validate = validate
        self.description = description
#        self.optional = optional
        self.dest = dest
        #self.default = default
        

    def __repr__(self):
        return 'Product(type=%r, dest=%r)' % (self.type, self.dest)


class Optional(object):
    def __init__(self, product_type):

        if inspect.isclass(product_type):
            product_type = product_type()

        if isinstance(product_type, DataProduct):
            self.type = product_type
        else:
            raise TypeError('product_type must be of class DataProduct')

class BaseRecipeResult(object):
    def __new__(cls, *args, **kwds):
        return super(BaseRecipeResult, cls).__new__(cls)

    def __init__(self, *args, **kwds):
        super(BaseRecipeResult, self).__init__()
    
    def suggest_store(self, *args, **kwds):
        pass

class ErrorRecipeResult(BaseRecipeResult):
    def __init__(self, errortype, message, traceback):
        self.errortype = errortype
        self.message = message
        self.traceback = traceback

    def __repr__(self):
        sclass = type(self).__name__
        return "%s(errortype=%r, message='%s')" % (sclass, 
            self.errortype, self.message)

class RecipeResultType(MapStoreType):
    '''Metaclass for RecipeResult.'''
    @classmethod
    def exclude(cls, name, value):
        return isinstance(value, Product)
    
    @classmethod
    def store(cls, name, value):
        if value.dest is None:
            value.dest = name
        nname = value.dest
        return nname, value

class RecipeResultAutoQCType(RecipeResultType):
    '''Metaclass for RecipeResult with added QC'''
    def __new__(cls, classname, parents, attributes):
        if 'qc' not in attributes:
            attributes['qc'] = Product(QualityControlProduct)
        return super(RecipeResultAutoQCType, cls).__new__(cls, classname, parents, attributes)


class RecipeResult(BaseRecipeResult):
    __metaclass__ = RecipeResultType

    def __new__(cls, *args, **kwds):
        self = super(RecipeResult, cls).__new__(cls)
        for key, prod in cls.iteritems():
            if key in kwds:
                val = kwds[key]
                nval = prod.type.store(val)
            elif prod.type.default:
                val = prod.type.default
                nval = prod.type.store(val)
            elif prod.optional:
                nval = None
            else:
                raise ValueError('required DataProduct %r not defined' %
                                 prod.type.__class__.__name__)
                
            setattr(self, key, nval)                
        return self

    def __init__(self, *args, **kwds):
        super(RecipeResult, self).__init__(self, *args, **kwds)

    def __repr__(self):
        sclass = type(self).__name__
        full = []
        for key, val in self.__class__.iteritems():
            full.append('%s=%r' % (key, val))
        return '%s(%s)' % (sclass, ', '.join(full))

    def suggest_store(self, **kwds):
        for k in kwds:
            mm = getattr(self, k)
            self.__class__[k].type.suggest(mm, kwds[k])

    def validate(self):
        '''Validate myself.'''

        # By default, validate each value
        for key, req in self.__class__.items():
            val = getattr(self, key)
            req.type.validate(val)

class RecipeResultAutoQC(RecipeResult):
    '''RecipeResult with an automatic QC member.'''
    __metaclass__ = RecipeResultAutoQCType

def transmit(result):
    if not isinstance(result, BaseRecipeResult):
        raise TypeError('result must be a RecipeResult')
    if isinstance(result, RecipeResult):
        pass # transmit as valid'
    elif isinstance(result, ErrorRecipeResult):
        res = {'error': {'type': result.errortype,
                         'message': result.message,
                         'traceback': result.traceback}
                         }
        return res
    else:
        raise TypeError('Unknown subclass of RecipeResult')

class define_result(object):
    def __init__(self, resultClass):
        if not issubclass(resultClass, BaseRecipeResult):
            raise TypeError('%r does not derive from BaseRecipeResult' % resultClass)

        self.klass = resultClass

    def __call__(self, klass):
        klass.RecipeResult = self.klass
        return klass

