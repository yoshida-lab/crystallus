import pytest

from crystallus import WyckoffCfgGenerator


def test_should_raise_value_err():
    with pytest.raises(ValueError, match="no configurations for generation"):
        WyckoffCfgGenerator()


def test_wyckoff_gen_init():
    wy = WyckoffCfgGenerator(Ca=2, C=2, O=6)
    assert wy.composition == dict(Ca=2, C=2, O=6)
    assert wy.priority == None

    wy = WyckoffCfgGenerator(Ca=2, C=2, O=6, priority={167: {'e': 0}})
    assert wy.composition == dict(Ca=2, C=2, O=6)
    assert wy.priority == {167: {'e': 0}}


def test_wyckoff_gen_without_priority():
    wy = WyckoffCfgGenerator(Ca=2, C=2, O=6)
    cfg = wy.gen_one(167)
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

    cfgs = wy.gen_many(1000, 167)
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

    cfgs = wy.gen_many(1000, 167, 166)
    assert isinstance(cfgs, dict)
    assert len(cfgs) == 2
    assert 166 in cfgs
    assert 167 in cfgs


def test_wyckoff_gen_with_priority():
    wy = WyckoffCfgGenerator(Ca=2, C=2, O=6, priority={167: {'e': 0}})
    cfg = wy.gen_one(167)
    assert cfg in [{'Ca': ['b'], 'C': ['a'], 'O': ['d']}, {'Ca': ['a'], 'C': ['b'], 'O': ['d']}]
    cfgs = wy.gen_many(1000, 167)
    assert all([
        cfg in cfgs for cfg in [{
            'Ca': ['b'],
            'C': ['a'],
            'O': ['d']
        }, {
            'Ca': ['a'],
            'C': ['b'],
            'O': ['d']
        }]
    ])
