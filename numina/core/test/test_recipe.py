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


'''Unit test for RecipeBase.'''

from six import with_metaclass

from ..metarecipes import RecipeType
from ..recipes import BaseRecipeAutoQC
from ..recipeinout import RecipeRequirements, RecipeResult
from ..recipeinout import RecipeResultAutoQC
from ..requirements import ObservationResultRequirement
from ..dataholders import Product
from ..products import QualityControlProduct


def test_metaclass_empty_base():

    class TestRecipe(with_metaclass(RecipeType, object)):
        pass

    assert hasattr(TestRecipe, 'RecipeRequirements')

    assert hasattr(TestRecipe, 'RecipeResult')

    assert hasattr(TestRecipe, 'Requirements')

    assert hasattr(TestRecipe, 'Result')

    assert issubclass(TestRecipe.RecipeRequirements, RecipeRequirements)

    assert issubclass(TestRecipe.RecipeResult, RecipeResult)

    assert TestRecipe.RecipeRequirements.__name__ == 'RecipeRequirements'

    assert TestRecipe.RecipeResult.__name__ == 'RecipeResult'


def test_metaclass():

    class TestRecipe(with_metaclass(RecipeType, object)):
        pass   

        obsresult = ObservationResultRequirement()
        someresult = Product(int, 'Some integer')

    assert hasattr(TestRecipe, 'RecipeRequirements')
    assert hasattr(TestRecipe, 'Requirements')

    assert hasattr(TestRecipe, 'RecipeResult')
    assert hasattr(TestRecipe, 'Result')

    assert issubclass(TestRecipe.Requirements, RecipeRequirements)

    assert issubclass(TestRecipe.Result, RecipeResult)

    assert TestRecipe.RecipeRequirements.__name__ == 'TestRecipeRequirements'

    assert TestRecipe.RecipeResult.__name__ == 'TestRecipeResult'


def test_recipe_with_autoqc():

    class TestRecipe(BaseRecipeAutoQC):
        obsresult = ObservationResultRequirement()
        someresult = Product(int, 'Some integer')

    assert hasattr(TestRecipe, 'RecipeRequirements')

    assert hasattr(TestRecipe, 'RecipeResult')

    assert hasattr(TestRecipe, 'Requirements')

    assert hasattr(TestRecipe, 'Result')

    assert issubclass(TestRecipe.RecipeRequirements, RecipeRequirements)

    assert issubclass(TestRecipe.RecipeResult, RecipeResultAutoQC)

    assert TestRecipe.RecipeRequirements.__name__ == 'TestRecipeRequirements'

    assert TestRecipe.RecipeResult.__name__ == 'TestRecipeResult'

    assert 'qc' in TestRecipe.RecipeResult

    for prod in TestRecipe.RecipeResult.values():
        assert isinstance(prod, Product)

    qc = TestRecipe.RecipeResult['qc']

    assert isinstance(qc.type, QualityControlProduct)
