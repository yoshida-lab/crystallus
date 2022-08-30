from crystallus import WyckoffCfgGenerator


def test_wyckoff_gen_init():
    wy = WyckoffCfgGenerator(dict(Ca=2, C=2, O=6))
    assert wy.composition == dict(Ca=2, C=2, O=6)
    assert wy.priority == None

    wy = WyckoffCfgGenerator(dict(Ca=2, C=2, O=6), priority={167: {'e': 0}})
    assert wy.composition == dict(Ca=2, C=2, O=6)
    assert wy.priority == {167: {'e': 0}}


def test_wyckoff_gen_without_priority():
    wy = WyckoffCfgGenerator(dict(Ca=2, C=2, O=6), verbose=False)

    cfg = wy.gen_one(spacegroup_num=167)
    assert cfg in [{
        'Ca': ['b'],
        'C': ['a'],
        'O': ['e']
    }, {
        'Ca': ['b'],
        'C': ['a'],
        'O': ['d']
    }, {
        'Ca': ['a'],
        'C': ['b'],
        'O': ['e']
    }, {
        'Ca': ['a'],
        'C': ['b'],
        'O': ['d']
    }]

    cfgs = wy.gen_many(1000, spacegroup_num=167)
    assert all([
        cfg in cfgs for cfg in [{
            'Ca': ['b'],
            'C': ['a'],
            'O': ['e']
        }, {
            'Ca': ['b'],
            'C': ['a'],
            'O': ['d']
        }, {
            'Ca': ['a'],
            'C': ['b'],
            'O': ['e']
        }, {
            'Ca': ['a'],
            'C': ['b'],
            'O': ['d']
        }]
    ])

    cfgs = wy.gen_many(1000, spacegroup_num=(167, 166))
    assert isinstance(cfgs, dict)
    assert len(cfgs) == 2
    assert 166 in cfgs
    assert 167 in cfgs


def test_wyckoff_gen_with_priority():
    wy = WyckoffCfgGenerator(dict(Ca=2, C=2, O=6), priority={167: {'e': 0, 'd': 2, 'b': 2, 'a': 2, 'c': 2, 'f': 2}})
    cfg = wy.gen_one(spacegroup_num=167)
    assert cfg in [{'Ca': ['b'], 'C': ['a'], 'O': ['d']}, {'Ca': ['a'], 'C': ['b'], 'O': ['d']}]
    cfgs = wy.gen_many(1000, spacegroup_num=167)
    assert all([cfg in cfgs for cfg in [{'Ca': ['b'], 'C': ['a'], 'O': ['d']}, {'Ca': ['a'], 'C': ['b'], 'O': ['d']}]])
