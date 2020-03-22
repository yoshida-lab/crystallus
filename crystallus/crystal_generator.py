# Copyright 2020 TsumiNa. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

from .crystallus import CrystalGenerator as _CG
from ._util import dangerwrap
from typing import Tuple, Dict

__all__ = ["CrystalGenerator"]


class CrystalGenerator(object):

    def __init__(
        self,
        spacegroup_num: int,
        estimated_volume: float,
        estimated_variance: float,
        *,
        min_distance_tolerance: float = 0.15,
        angle_range: Tuple[float, float] = (30., 150.),
        angle_tolerance: float = 20.,
        max_recurrent: int = 5_000,
        n_jobs: int = -1,
    ):
        """A generator to generate crystal structures.

        Parameters
        ----------
        spacegroup_num : int
            Specify the spacegroup.
        estimated_volume : float
            The estimated volume of primitive cell.
        estimated_variance : float
            The estimated variance of volume prediction.
        min_distance_tolerance : float, optional
            The tolerance of atomic distances when distance checking, by default 0.15
        angle_range : Tuple[float, float], optional
            The range of angles when lattice generation, by default (30., 150.)
        angle_tolerance : float, optional
            The Tolerance of minimum angles when lattice generation, by default 20.
        max_recurrent : int, optional
            Max recurrent until generate a reasonable structure, by default 5_000
        n_jobs : int, optional
            Number of cpu cores when parallel calculation, by default -1
        """
        self._cg = _CG(spacegroup_num=spacegroup_num,
                       estimated_volume=estimated_volume,
                       estimated_variance=estimated_variance,
                       min_distance_tolerance=min_distance_tolerance,
                       angle_range=angle_range,
                       angle_tolerance=angle_tolerance,
                       max_recurrent=max_recurrent,
                       n_jobs=n_jobs)
        self._estimated_volume = estimated_volume
        self._estimated_variance = estimated_variance
        self._min_distance_tolerance = min_distance_tolerance
        self._angle_range = angle_range
        self._angle_tolerance = angle_tolerance
        self._max_recurrent = max_recurrent
        self._spacegroup_num = spacegroup_num

    @property
    def estimated_volume(self):
        return self._estimated_volume

    @property
    def estimated_variance(self):
        return self._estimated_variance

    @property
    def min_distance_tolerance(self):
        return self._min_distance_tolerance

    @property
    def angle_range(self):
        return self._angle_range

    @property
    def angle_tolerance(self):
        return self._angle_tolerance

    @property
    def max_recurrent(self):
        return self._max_recurrent

    @property
    def spacegroup_num(self):
        return self._spacegroup_num

    @property
    def n_jobs(self):
        return self._cg.n_jobs

    @n_jobs.setter
    def n_jobs(self, n):
        self._cg.n_jobs = n

    def gen_one(self, **cfg: Dict[str, Tuple[str]]):
        """Generate one structure with given configration set.

        Parameters
        ----------
        **cfg: Dict[str, Tuple[str]]
            Configuration set with format like: {"Li": ("a", "c"), "O": ("i",)},
            where "Li" is an available element symbol and ("a", "c") is a tuple
            contains coresponding Wyckoff letters.

        Returns
        -------
        Dict
            Structure information contains ``spacegroup_mun: int``,
            ``volume: float``, ``lattice: list``, ``wyckoff_letters: list``,
            and ``coords: list``.
        """
        return self._cg.gen_one(**cfg)

    def gen_many(self,
                 size: int,
                 *cfgs: Dict[str, Tuple[str]],
                 iterative: bool = False,
                 **cfg: Dict[str, Tuple[str]]):
        """Generate structures with given configration set(s) and size.

        Parameters
        ----------
        *cfgs: Dict[str, Tuple[str]]
            A tuple with configuration set ``cfg``s.
        size: int
            Generation size for each configration set. By default, 1
        iterative: bool
            Running like a iterator. Instead of return results until all done,
            Results will be returned when generation of each configration has done.
            This only works when length of ``cfgs`` greater than 0. By default, False
        **cfg: Dict[str, Tuple[str]]
            Configuration set with format like: {"Li": ("a", "c"), "O": ("i",)},
            where "Li" is an available element symbol and ("a", "c") is a tuple
            which contains coresponding Wyckoff letters. For convenience, dict will
            be sorted by keys.

        Returns
        -------
        Dict
            Structure information contains ``spacegroup_mun: int``,
            ``volume: float``, ``lattice: list``, ``wyckoff_letters: list``,
            and ``coords: list``.

        Yields
        ------
        Tuple[Dict]
            Structure information contains ``spacegroup_mun: int``,
            ``volume: float``, ``lattice: list``, ``wyckoff_letters: list``,
            and ``coords: list``.
        """
        assert size >= 1, 'size must be greater than 1'
        if len(cfgs) != 0 and len(cfg) != 0:
            raise ValueError('`cfgs` and `cfg are exclusive')

        if len(cfg) != 0:
            return dangerwrap(self._cg.gen_many, size, cfg)
        if len(cfgs) != 0:
            if iterative:
                for cfg in cfgs:
                    yield dangerwrap(self._cg.gen_many, size, cfg)
            else:
                return dangerwrap(self._cg.gen_many, size, *cfgs)

        raise ValueError('need configration dicts as `cfgs` or dict as `cfg')

    def __repr__(self):
        return f"CrystalGenerator(\
            \n    spacegroup_num={self.spacegroup_num},\
            \n    estimated_volume={self.estimated_volume},\
            \n    estimated_variance={self.estimated_variance},\
            \n     min_distance_tolerance={self.min_distance_tolerance},\
            \n    angle_range={self.angle_range},\
            \n    angle_tolerance={self.angle_tolerance},\
            \n    max_recurrent={self.max_recurrent},\
            \n    n_jobs={self.n_jobs}\
            \n)"
