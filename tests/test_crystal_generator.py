# Copyright 2021 TsumiNa
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pytest
import numpy as np
import numpy.testing as npt

from crystallus import CrystalGenerator


def test_crystal_gen_one():
    cg = CrystalGenerator(207, 1000, 20)
    structure = cg.gen_one(C=('a', 'b'), O=('d'))

    assert set(structure.keys()) == {
        "lattice",
        "volume",
        "spacegroup_num",
        "species",
        "wyckoff_letters",
        "coords",
    }
    assert 940 <= structure['volume'] <= 1060
    l = pow(structure['volume'], 1 / 3)
    assert np.allclose(structure['lattice'], [[l, 0, 0], [0, l, 0], [0, 0, l]])
    assert structure['spacegroup_num'] == 207
    assert structure['species'] == ['C', 'C', 'O', 'O', 'O']
    assert structure['wyckoff_letters'] == ['a', 'b', 'd', 'd', 'd']
    assert structure['coords'] == [[0., 0., 0.], [1 / 2, 1 / 2, 1 / 2], [1 / 2, 0., 0.],
                                   [0., 1 / 2, 0.], [0., 0., 1 / 2]]


def test_crystal_gen_one_with_template():
    """
    space group 167
    =========================================================================
    Multiplicity | Wyckoff letter |      Coordinates
    -------------------------------------------------------------------------
            6       |         e      | (x,0,1/4) (0,x,1/4) (-x,-x,1/4)
                    |                | (-x,0,3/4) (0,-x,3/4) (x,x,3/4)
    -------------------------------------------------------------------------
            6       |         d      | (1/2,0,0) (0,1/2,0) (1/2,1/2,0)
                    |                | (0,1/2,1/2) (1/2,0,1/2) (1/2,1/2,1/2)
    -------------------------------------------------------------------------
            4       |         c      | (0,0,z) (0,0,-z+1/2) (0,0,-z) (0,0,z+1/2)
    -------------------------------------------------------------------------
            2       |         b      | (0,0,0) (0,0,1/2)
    -------------------------------------------------------------------------
            2       |         a      | (0,0,1/4) (0,0,3/4)
    -------------------------------------------------------------------------
    """
    template = [('c', [0.2, 0., 0.0]), ('e', [0.4, 0., 0.0])]
    cg = CrystalGenerator(167, 1000, 10, empirical_coords=template, empirical_coords_variance=0)
    structure = cg.gen_one(C=('c'), O=('e'))

    assert structure['wyckoff_letters'] == ['c'] * 4 + ['e'] * 6
    npt.assert_almost_equal(
        structure['coords'],
        [
            [0.2, 0.2, 0.2],  # Li
            [0.3, 0.3, 0.3],  # Li
            [0.8, 0.8, 0.8],  # Li
            [0.7, 0.7, 0.7],  # Li
            [0.4, 0.1, 0.25],  # P
            [0.25, 0.4, 0.1],  # P
            [0.1, 0.25, 0.4],  # P
            [0.6, 0.9, 0.75],  # P
            [0.75, 0.6, 0.9],  # P
            [0.9, 0.75, 0.6],  # P
        ])


def test_crystal_gen_many_1():
    cg = CrystalGenerator(207, 1000, 20)
    structure = cg.gen_many(10)

    assert len(structure) == 0


def test_gen_many_2():
    cg = CrystalGenerator(207, 1000, 20)
    with pytest.raises(ValueError, match="`max_attempts` can not be smaller than `expect_size`"):
        cg.gen_many(10, {'C': ('a', 'b'), 'O': ('d',)}, max_attempts=2)


def test_crystal_gen_many_3():
    cg = CrystalGenerator(207, 1000, 20)
    structure = cg.gen_many(10, {'C': ('a', 'b'), 'O': ('d',)})

    assert len(structure) == 10
    assert structure[0]['wyckoff_letters'] == structure[2]['wyckoff_letters']
    assert structure[0]['species'] == structure[2]['species']
    assert structure[0]['spacegroup_num'] == structure[2]['spacegroup_num']
    assert structure[0]['coords'] == structure[2]['coords']


def test_crystal_gen_many_4():
    cg = CrystalGenerator(207, 1000, 20)
    cfgs = (
        {
            'C': ('a', 'b'),
            'O': ('d',)
        },
        {
            'O': ('d',),
            'C': ('b', 'a')
        },
    )
    structure = cg.gen_many(5, *cfgs)

    assert len(structure) == 10
    assert structure[0]['wyckoff_letters'] == ['a', 'b', 'd', 'd', 'd']
    assert structure[5]['wyckoff_letters'] == ['b', 'a', 'd', 'd', 'd']
    assert structure[0]['wyckoff_letters'] == structure[4]['wyckoff_letters']
    assert structure[5]['wyckoff_letters'] == structure[9]['wyckoff_letters']


def test_crystal_gen_many_5():
    cg = CrystalGenerator(33, 1168, 15)
    comp = {
        'Ag': ['a', 'a', 'a', 'a', 'a', 'a', 'a', 'a'],
        'Ge': ['a'],
        'S': ['a', 'a', 'a', 'a', 'a', 'a']
    }

    # distance error
    structure = cg.gen_many(100, comp)
    assert len(structure) == 0

    # no check
    structure = cg.gen_many(10, comp, check_distance=False)
    assert len(structure) == 10

    # make condition losser
    structure = cg.gen_many(1000, comp, distance_scale_factor=0.5)
    assert len(structure) > 0
