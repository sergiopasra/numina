#
# Copyright 2014-2017 Universidad Complutense de Madrid
#
# This file is part of Numina
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""DAL base classes"""


from .daliface import DALInterface
from numina.exceptions import NoResultFound


class AbsDAL(DALInterface):
    pass


class AbsDrpDAL(DALInterface):
    def __init__(self, drps, *args, **kwargs):
        super(AbsDrpDAL, self).__init__()
        self.drps = drps

    def search_instrument_configuration_from_ob(self, ob):
        ins = ob.instrument
        name = ob.configuration
        return self.search_instrument_configuration(ins, name)

    def search_instrument_configuration(self, ins, name):

        drp = self.drps.query_by_name(ins)

        try:
            this_configuration = drp.configurations[name]
        except KeyError:
            raise KeyError('Instrument configuration "{}" missing'.format(name))

        return this_configuration

    def search_recipe(self, ins, mode, pipeline):

        drp = self.drps.query_by_name(ins)

        if drp is None:
            raise NoResultFound('DRP not found')

        try:
            this_pipeline = drp.pipelines[pipeline]
        except KeyError:
            raise NoResultFound('pipeline not found')

        try:
            recipe = this_pipeline.get_recipe_object(mode)
            return recipe
        except KeyError:
            raise NoResultFound('mode not found')

    def search_recipe_fqn(self, ins, mode, pipename):

        drp = self.drps.query_by_name(ins)

        this_pipeline = drp.pipelines[pipename]
        recipes = this_pipeline.recipes
        recipe_fqn = recipes[mode]
        return recipe_fqn

    def search_recipe_from_ob(self, ob, pipeline='default'):
        instrument = ob.instrument
        mode = ob.mode
        pipeline = ob.pipeline
        return self.search_recipe(instrument, mode, pipeline)
