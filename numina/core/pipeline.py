#
# Copyright 2011-2014 Universidad Complutense de Madrid
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

'''DRP loader.'''

import logging

import pkg_resources


_logger = logging.getLogger('numina')


class Pipeline(object):
    '''Base class for pipelines.'''
    def __init__(self, name, recipes, version=1):
        self.name = name
        self.recipes = recipes
        self.version = version

    def get_recipe(self, mode):
        return self.recipes[mode]


class InstrumentConfiguration(object):
    '''Configuration of an Instrument.'''
    def __init__(self, values):
        self.values = values


class Instrument(object):
    '''Description of an Instrument.'''
    def __init__(self, name, configurations, modes, pipelines):
        self.name = name
        self.configurations = configurations
        self.modes = modes
        self.pipelines = pipelines


class ObservingMode(object):
    '''Observing modes of an Instrument.'''
    def __init__(self):
        self.name = ''
        self.uuid = ''
        self.key = ''
        self.url = ''
        self.instrument = ''
        self.summary = ''
        self.description = ''
        self.recipe = ''
        self.recipe_class = None
        self.status = ''
        self.date = ''
        self.reference = ''


class LoadableDRP(object):
    '''Container for the loaded DRP.'''
    def __init__(self, instruments):
        self.instruments = instruments


def init_drp_system():
    '''Load all available DRPs in 'numina.pipeline' entry_point.'''

    drp = {}

    for entry in pkg_resources.iter_entry_points(group='numina.pipeline.1'):
        drp_loader = entry.load()
        mod = drp_loader()
        if mod:
            drp.update(mod.instruments)
        else:
            _logger.warning('Module %s does not contain a valid DRP', mod)

    return drp


def init_backends(backend='default'):
    '''Load storage modes'.'''

    backends = []

    for entry in pkg_resources.iter_entry_points(group='numina.storage.1'):
        store = entry.load()
        backends.append((store, True))

    return backends
