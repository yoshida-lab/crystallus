import pytest
import numpy as np

from crystallus import CrystalGenerator


def test_gen_one():
    cg = CrystalGenerator(207, 1000, 20)
    structure = cg.gen_one(C=('a', 'b'), O=('d'))

    assert set(structure.keys()) == {
        "lattice", "volume", "spacegroup_num", "species", "wyckoff_letters", "coords",
        "attempts_until_done"
    }
    assert 940 <= structure['volume'] <= 1060
    l = pow(structure['volume'], 1 / 3)
    assert np.allclose(structure['lattice'], [[l, 0, 0], [0, l, 0], [0, 0, l]])
    assert structure['spacegroup_num'] == 207
    assert structure['species'] == ['C', 'C', 'O', 'O', 'O']
    assert structure['wyckoff_letters'] == ['a', 'b', 'd', 'd', 'd']
    assert structure['coords'] == [[0., 0., 0.], [1 / 2, 1 / 2, 1 / 2], [1 / 2, 0., 0.],
                                   [0., 1 / 2, 0.], [0., 0., 1 / 2]]


def test_gen_many_1():
    cg = CrystalGenerator(207, 1000, 20)
    structure = cg.gen_many(10, {'C': ('a', 'b'), 'O': ('d',)})

    assert len(structure) == 10
    assert structure[0]['wyckoff_letters'] == structure[2]['wyckoff_letters']
    assert structure[0]['species'] == structure[2]['species']
    assert structure[0]['spacegroup_num'] == structure[2]['spacegroup_num']
    assert structure[0]['coords'] == structure[2]['coords']


def test_gen_many_2():
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


def test_gen_many_3():
    cg = CrystalGenerator(33, 1168, 15)
    comp = {
        'Ag': ['a', 'a', 'a', 'a', 'a', 'a', 'a', 'a'],
        'Ge': ['a'],
        'S': ['a', 'a', 'a', 'a', 'a', 'a']
    }

    # distance error
    structure = cg.gen_many(10, comp)
    assert len(structure) == 0

    # no check
    structure = cg.gen_many(10, comp, check_distance=False)
    assert len(structure) == 10

    # make condition losser
    structure = cg.gen_many(10, comp, atomic_distance_tolerance=0.5)
    assert len(structure) > 0
