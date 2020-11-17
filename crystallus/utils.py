import re
from copy import deepcopy
from typing import List, Tuple, Union, Callable

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
    """Convert fraction coordinates into Wyckoff position formate. """

    patten = re.compile(r'(?<=\)),\s*')

    def __init__(self, spacegroup_num: int):
        wys = SpaceGroupDB.get(SpaceGroupDB.spacegroup_num == spacegroup_num).wyckoffs
        self.particle = _Particle()
        self.wyckoff_pos = {wy.letter: self.patten.split(wy.positions)[0][1:-1] for wy in wys}

    def _inner(self, wy_letter, coord):
        b = np.asarray(coord)
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

    def __call__(self,
                 wyckoff_letters: Union[str, pd.Series, List[str]],
                 coords: Union[str, pd.Series, List[Tuple[float, float, float]]],
                 elements: Union[str, pd.Series, List[str], None] = None,
                 *,
                 data: pd.DataFrame = None):
        if data is not None:
            if not isinstance(wyckoff_letters, str) or not isinstance(coords, str):
                raise ValueError(
                    '`wyckoff_letters` and `coords` must be the column name when `data` is set')

            if elements is not None and isinstance(elements, str):
                wy_and_coord = [a for _, a in data[[wyckoff_letters, coords, elements]].iterrows()]
            else:
                wy_and_coord = [a for _, a in data[[wyckoff_letters, coords]].iterrows()]
        else:
            if isinstance(wyckoff_letters, str) or isinstance(coords, str) or isinstance(
                    elements, str):
                raise ValueError(
                    'found `wyckoff_letters`, `coords`, and `elements` as column name but `data` is not given'
                )

            if elements is None:
                if not len(wyckoff_letters) == len(coords):
                    raise ValueError('`wyckoff_letters`and `coords` have different lengths')
                wy_and_coord = list(zip(wyckoff_letters, coords))
            else:
                if not len(wyckoff_letters) == len(coords) == len(elements):
                    raise ValueError(
                        '`wyckoff_letters`, `coords`, and `elements` have different lengths')
                wy_and_coord = list(zip(wyckoff_letters, coords, elements))

        if len(wy_and_coord[0]) == 2:
            return [(wy, self._inner(wy, b)) for wy, b in wy_and_coord]
        if len(wy_and_coord[0]) == 3:
            return [(f'{elem}:{wy}', self._inner(wy, b)) for wy, b, elem in wy_and_coord]
        raise ValueError(
            '`wy_and_coord` must be a list of (wyckoff_letter, coord) or (wyckoff_letter, coord, element)'
        )


def build_structure(structure_data: dict):
    """Build a pymatgen Structure object from the genrated structure data.

    Parameters
    ----------
    structure_data:
        A generate structure data.

    Returns
    -------
    dict
        A dictionary contains a structure and other information about the structure.
    """
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


def get_equivalent_coords(structure: Structure, *, mapper: Callable[[str, str, int], str] = None):
    """Extract the equivalent coordinates from the given structure.

    Parameters
    ----------
    structure:
        A pymatgen structure object.
    mapper:
        Specify how to replace the elements. optional.
        ``mapper`` has signature ``[element, wyckoff_letter, multiplicity] -> target_element``
        If this parameter is given, will map the element in the structure to the corresponding one.
        For example, replace the `Ca` in `CaCO2` with `Mg`.

    Returns
    -------
    DataFrame
        A dataframe contains all equivalent coordinates and their Wyckoff position letters.
    """
    struct = SpacegroupAnalyzer(structure).get_symmetrized_structure()

    def _inner(i, sites, mapper=None):
        site = sites[0]
        wy_symbol = struct.wyckoff_symbols[i]
        row = {'element': site.species_string}
        row['spacegroup_num'] = struct.get_space_group_info()[1]
        row['multiplicity'] = int(wy_symbol[:-1])
        row['wyckoff_letter'] = wy_symbol[-1]
        if mapper is not None:
            row['target_element'] = mapper(site.species_string, row['wyckoff_letter'],
                                           row['multiplicity'])
        row['coordinate'] = list(site.frac_coords)

        return row

    return pd.DataFrame(
        [_inner(i, sites, mapper=mapper) for i, sites in enumerate(struct.equivalent_sites)])


def structure_dissimilarity(anchor_structure: Structure,
                            other_structures: List[Structure],
                            *,
                            verbose: int = 1,
                            n_jobs: int = 1):
    """Calculate dissimilarity between anchor and other structures.

    Parameters
    ----------
    anchor_structure:
        Anchor structure
    other_structures:
        Structures will be used to calculate the dissimilarity against the anchor structure.
    verbose:
        Verbose output when performing parallel calculation, by default 1
    n_jobs:
        Specify the number of cores for parallel calculation, by default 1

    Returns
    -------
    list
        List of dissimilarities.
    """
    # Calculate structure fingerprints.
    ssf = SiteStatsFingerprint(CrystalNNFingerprint.from_preset('ops',
                                                                distance_cutoffs=None,
                                                                x_diff_weight=0),
                               stats=('mean', 'std_dev', 'minimum', 'maximum'))
    v_anchor = np.array(ssf.featurize(anchor_structure))
    tmp = Parallel(n_jobs=n_jobs,
                   verbose=verbose)(delayed(ssf.featurize)(s) for s in other_structures)
    return [np.linalg.norm(np.array(s) - v_anchor) for s in tmp]
