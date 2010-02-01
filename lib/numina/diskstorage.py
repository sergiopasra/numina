#
# Copyright 2008-2010 Sergio Pascual
# 
# This file is part of PyEmir
# 
# PyEmir is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# PyEmir is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PyEmir.  If not, see <http://www.gnu.org/licenses/>.
# 

# $Id$

__version__ = "$Revision$"

import simplejson as json

from pyfits.NP_pyfits import HDUList

_internal_map = {}

def _store_fits(obj):
    where = val['primary'].header.get('filename')
    obj.writeto(where, clobber=True, output_verify='ignore')
    
_internal_map[HDUList] = _store_fits 

def _store_map(obj):
    filename = 'products.json'
    f = open(filename, 'w+') 
    try:
        json.dump(obj, f)
    finally:
        f.close()

_internal_map[type({})] = _store_map

_internal_map[type([])] = _store_map

def register(cls):
    def wrap(f):
        _internal_map[cls] = f
        def w_f(*args):
            return f(*args)
        
        return w_f
    
    return wrap

def store_to_disk(result, filename='products.json'):
    rep = {}
    for key, val in result.products.iteritems():
        if val.__class__ in _internal_map:
            store = _internal_map[val.__class__]
            store(val)
        else:
            rep[key] = val
    
    f = open(filename, 'w+') 
    try:
        json.dump(rep, f)
    finally:
        f.close()
