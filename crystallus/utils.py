import re
from copy import deepcopy

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from matminer.featurizers.site import CrystalNNFingerprint
from matminer.featurizers.structure import SiteStatsFingerprint
from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from crystallus import SpaceGroupDB

__call__ = [
    'WyckoffPositionConverter', 'build_structure', 'get_equivalent_coords',
    'structure_dissimilarity'
]


class _Coordinate:
    patten = re.compile(r'(?P<xyz>-?\d?[xyz])|(?P<cons_frac>\d\/\d?)|(?P<cons>\d)')

    @classmethod
    def _inner(cls, s):
        if '-' in s:
            coeff = -1
        else:
            coeff = 1
        if '2' in s:
            coeff *= 2

        return coeff

    def __call__(self, coordinates):
        const = 0
        x_coeff, y_coeff, z_coeff, const = 0, 0, 0, 0

        for e in self.patten.findall(coordinates):

            if e[0] != '':
                s = e[0].lower()
                if 'x' in s:
                    x_coeff = self._inner(s)
                    continue

                if 'y' in s:
                    y_coeff = self._inner(s)
                    continue

                if 'z' in s:
                    z_coeff = self._inner(s)
                    continue

            if e[1] != '':
                s = e[1].split('/')
                const = float(s[0]) / float(s[1])
                continue

            if e[2] != '':
                const = float(e[2])
                continue

        return x_coeff, y_coeff, z_coeff, const


class _Particle():
    patten = re.compile(r',\s*')

    def __init__(self):
        self.Coordinate = _Coordinate()

    def __call__(self, position):
        return [self.Coordinate(coord) for coord in self.patten.split(position)]


class WyckoffPositionConverter():
    patten = re.compile(r'(?<=\)),\s*')

    def __init__(self, spacegroup_num):
        wys = SpaceGroupDB.get(SpaceGroupDB.spacegroup_num == spacegroup_num).wyckoffs
        self.particle = _Particle()
        self.wyckoff_pos = {wy.letter: self.patten.split(wy.positions)[0][1:-1] for wy in wys}

    def __call__(self, wy_letter, b):
        b = np.asarray(b)
        a = np.array(self.particle(self.wyckoff_pos[wy_letter]))
        idx = []

        if np.count_nonzero(a[:, 0]):
            idx.append(0)
        if np.count_nonzero(a[:, 1]):
            idx.append(1)
        if np.count_nonzero(a[:, 2]):
            idx.append(2)
        b[idx] -= a[idx, -1]

        if len(idx) > 1:
            solves = np.linalg.solve(a[idx][:, idx], b[idx] - a[idx, -1])
            b[idx] = solves

        return b.tolist()


def build_structure(structure_data):
    s = deepcopy(structure_data)
    struct = Structure(lattice=s['lattice'], species=s['species'], coords=s['coords'])
    s['structure'] = struct.get_primitive_structure()
    s['volume'] = struct.volume
    s['composition'] = dict(struct.composition.as_dict())
    s['formula'] = struct.composition.formula
    s['reduced_formula'] = struct.composition.reduced_formula
    s['num_atoms'] = struct.composition.num_atoms
    del s['coords']
    del s['lattice']

    return s


def get_equivalent_coords(structure: Structure):
    struct = SpacegroupAnalyzer(structure).get_symmetrized_structure()

    def _inner(i, sites):
        site = sites[0]
        wy_symbol = struct.wyckoff_symbols[i]
        row = {'element': site.species_string}
        row['multiplicity'] = int(wy_symbol[:-1])
        row['wyckoff_letter'] = wy_symbol[-1]
        row['coordinate'] = [j for j in site.frac_coords]

        return row

    return pd.DataFrame([_inner(i, sites) for i, sites in enumerate(struct.equivalent_sites)])


def structure_dissimilarity(anchor_structure, other_structures, *, verbose=1, n_jobs=1):
    # Calculate structure fingerprints.
    ssf = SiteStatsFingerprint(CrystalNNFingerprint.from_preset('ops',
                                                                distance_cutoffs=None,
                                                                x_diff_weight=0),
                               stats=('mean', 'std_dev', 'minimum', 'maximum'))
    v_anchor = np.array(ssf.featurize(anchor_structure))
    tmp = Parallel(n_jobs=n_jobs,
                   verbose=verbose)(delayed(ssf.featurize)(s) for s in other_structures)
    return [np.linalg.norm(np.array(s) - v_anchor) for s in tmp]
