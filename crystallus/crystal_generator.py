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

from .crystallus import CrystalGenerator as _CG
from copy import deepcopy
from typing import Tuple, Dict, List, Union

__all__ = ["CrystalGenerator"]


class CrystalGenerator(object):

    def __init__(self,
                 spacegroup_num: int,
                 estimated_volume: float,
                 estimated_variance: float,
                 *,
                 angle_range: Tuple[float, float] = (30., 150.),
                 angle_tolerance: float = 20.,
                 empirical_coords: Union[List[Tuple[str, List[float]]], None] = None,
                 empirical_coords_variance: float = 0.01,
                 empirical_coords_sampling_rate: float = 1.,
                 empirical_coords_loose_sampling: bool = True,
                 max_attempts_number: int = 5_000,
                 n_jobs: int = -1,
                 verbose: bool = False):
        """A generator for possible crystal structure generation.

        Parameters
        ----------
        spacegroup_num : int
            Specify the spacegroup.
        estimated_volume : float
            The estimated volume of primitive cell. Unit is Å^3.
        estimated_variance : float
            The estimated variance of volume prediction. Unit is Å^3.
            ``estimated_volume`` and ``estimated_variance`` will be used to build
            a Gaussion distribution for the sampling of volume of primitive cell.
            We will use the abstract valuse if sampled volume is negative.
        angle_range : Tuple[float, float], optional
            The range of the degree of angles when lattice generation. by default (30., 150.)
        angle_tolerance : float, optional
            The Tolerance of minimum of the degree of angles when lattice generation, by default 20.
        empirical_coords:
            Empirical distributuion of atomic coordinates. The coordinates should be given in Wyckoff position
            format. For example: for Wyckoff letter `c` in space group 167, the corresponding positions are
            `(0,0,z) (0,0,-z+1/2) (0,0,-z) (0,0,z+1/2)`. So a fraction coordinate in lattice such as
            `[[0,0,0.3] [0,0,0.2] [0,0,0.7] [0,0,0.8]]` should be converted to `([0, 0, 0.3])`.
            The key of each coordinates pool has two mode, the one is using Wyckoff position letter, for example
            `[('a', [0, 0, 0.3]), ('b': [...])]`, the other is by element and Wyckoff position letter,
            for example `[('Ca:a': [0, 0, 0.3]), ('Ca:b': [...])], ('C:a': [...])]`. You can use parameter
            `empirical_coords_loose_sampling` to switch modes.
        empirical_coords_variance:
            The variance of empirical_coords. This parameter will be used to build a Gaussian distribution.
            The generator will sample values from the distribution as the perturbation of empirical coordinates.
        empirical_coords_sampling_rate:
            The sampling rate when sampling the empirical coordinates.
            A higher rate means instead of random generation,
            sampling from empirical coordinate distribution as more as possible when generating.
            Must be in 0 ~ 1, default 1.
        empirical_coords_loose_sampling:
            Specify the sampling mode. `True` for using Wyckoff position letter, `False` means using
            element and Wyckoff position letter. See parameter `empirical_coords` for details.
        max_attempts_number : int, optional
            Max recurrent until generate a reasonable lattice, by default is 5_000
        n_jobs : int, optional
            Number of cpu cores when parallel calculation, by default -1
        verbose: bool, optional
            Set to ``True`` to show more information.
        """
        # pyo3 can not convert PyList to PyTuple automatically
        if empirical_coords is not None and not isinstance(empirical_coords, tuple):
            empirical_coords = tuple(empirical_coords)

        self._estimated_volume = estimated_volume
        self._estimated_variance = estimated_variance
        self._angle_range = angle_range
        self._angle_tolerance = angle_tolerance
        self._max_attempts_number = max_attempts_number
        self._empirical_coords = empirical_coords
        self._empirical_coords_variance = empirical_coords_variance
        self._empirical_coords_sampling_rate = empirical_coords_sampling_rate
        self._empirical_coords_loose_sampling = empirical_coords_loose_sampling
        self._spacegroup_num = spacegroup_num

        self._cg = _CG(spacegroup_num=spacegroup_num,
                       estimated_volume=estimated_volume,
                       estimated_variance=estimated_variance,
                       angle_range=angle_range,
                       angle_tolerance=angle_tolerance,
                       empirical_coords=empirical_coords,
                       empirical_coords_variance=empirical_coords_variance,
                       empirical_coords_sampling_rate=empirical_coords_sampling_rate,
                       empirical_coords_loose_sampling=empirical_coords_loose_sampling,
                       max_attempts_number=max_attempts_number,
                       n_jobs=n_jobs,
                       verbose=verbose)

    @property
    def estimated_volume(self):
        return self._estimated_volume

    @property
    def estimated_variance(self):
        return self._estimated_variance

    @property
    def angle_range(self):
        return self._angle_range

    @property
    def angle_tolerance(self):
        return self._angle_tolerance

    @property
    def max_attempts_number(self):
        return self._max_attempts_number

    @property
    def empirical_coords(self):
        return deepcopy(self._empirical_coords)

    @property
    def empirical_coords_variance(self):
        return self._empirical_coords_variance

    @property
    def empirical_coords_sampling_rate(self):
        return self._empirical_coords_sampling_rate

    @property
    def empirical_coords_loose_sampling(self):
        return self._empirical_coords_loose_sampling

    @property
    def spacegroup_num(self):
        return self._spacegroup_num

    @property
    def n_jobs(self):
        return self._cg.n_jobs

    @n_jobs.setter
    def n_jobs(self, n):
        self._cg.n_jobs = n

    @property
    def verbose(self):
        return self._cg.verbose

    @verbose.setter
    def verbose(self, n):
        self._cg.verbose = n

    def gen_one(self,
                *,
                check_distance: bool = True,
                distance_scale_factor: float = 0.1,
                **cfg: Dict[str, Tuple[str]]):
        """Try to generate a legal crystal structure with given configuration set.

        Parameters
        ----------
        check_distance: bool, optional
            Whether the atomic distance should be checked. default ``True``
        distance_scale_factor : float, optional
            Scale factor to determine the tolerance of atomic distances when distance checking. Unit is Å,
            When ``check_distance`` is ``True``, Any structure has
            all_atomic_distance < (A_atom_covalent_radius + B_atom_covalent_radius) * (1 - distance_scale_factor) will be rejected,
            by default 0.1
        **cfg: Dict[str, Tuple[str]]
            Wyckoff Configuration set, which is a dict with format like:
            {"Li": ["a", "c"], "O": ["i"]}. Here, the "Li" is an available element
            symbol and ["a", "c"] is a list which contains coresponding Wyckoff
            letters. For convenience, dict will be sorted by keys.

        Returns
        -------
        Dict
            Structure information contains ``spacegroup_mun: int``,
            ``volume: float``, ``lattice: list``, ``wyckoff_letters: list``,
            and ``coords: list``.
        """
        return self._cg.gen_one(check_distance, distance_scale_factor, **cfg)

    def gen_many(
        self,
        expect_size: int,
        *cfgs: Dict[str, Tuple[str]],
        max_attempts: Union[int, None] = None,
        check_distance: bool = True,
        distance_scale_factor: float = 0.1,
    ) -> List[Dict]:
        """Try to generate legal crystal structures with given configuration set(s).

        Parameters
        ----------
        expect_size: int
            The expectation of the total amount of generated structures based on one Wyckoff.
            Whatever one generated structure is legal or not, **one attempt** will be consumed. 
            Please noted that the result could be empty when no structures matched the atomic distance conditions.
            When the number of generated structures are not fit your expectation too far away,
            try to give the parameter ``max_attempts`` a higher value..
        max_attempts: Union[int, None], optional
            Specify the max number of attempts in structure generation.
            When the number of generated structures is small than ``expect_size``, new rounds of structure generation will be performed.
            The generation will stop until the number of generated structures is more than ``expect_size, `` or the total attempts reach the ``max_attempts``.
            Default ``None``, means ``max_attempts`` equal to parameter ``expect_size``.
        check_distance: bool, optional
            Whether the atomic distance should be checked. default ``True``
        distance_scale_factor : float, optional
            Scale factor to determine the tolerance of atomic distances when distance checking. Unit is Å,
            When ``check_distance`` is ``True``, Any structure has
            all_atomic_distance < (A_atom_covalent_radius + B_atom_covalent_radius) * (1 - distance_scale_factor) will be rejected,
            by default 0.1
        *cfgs: Dict[str, Tuple[str]]
            A tuple with Wyckoff configuration set(s).
            Wyckoff Configuration set is a dict with format like: {"Li": ["a", "c"], "O": ["i"]}.
            Here, the "Li" is an available element symbol and ["a", "c"] is a list
            which contains coresponding Wyckoff letters. For convenience, dict will
            be sorted by keys..

        Returns
        -------
        Dict
            Structure information contains ``spacegroup_mun: int``,
            ``volume: float``, ``lattice: list``, ``wyckoff_letters: list``,
            and ``coords: list``.
        """
        assert expect_size >= 1, 'attempts number must be greater than 1'

        # if len(cfgs) > 0:
        return self._cg.gen_many(
            expect_size,
            max_attempts,
            check_distance,
            distance_scale_factor,
            *cfgs,
        )
        # return []

    def gen_many_iter(
        self,
        expect_size: int,
        *cfgs: Dict[str, Tuple[str]],
        max_attempts: Union[int, None] = None,
        check_distance: bool = True,
        distance_scale_factor: float = 0.1,
    ):
        """Try to generate legal crystal structures with given configuration set(s), iteratively.

        Parameters
        ----------
        expect_size: int
            The expectation of the total amount of generated structures based on one Wyckoff.
            Whatever one generated structure is legal or not, **one attempt** will be consumed. 
            Please noted that the result could be empty when no structures matched the atomic distance conditions.
            When the number of generated structures are not fit your expectation too far away,
            try to give the parameter ``max_attempts`` a higher value..
        max_attempts: Union[int, None], optional
            Specify the max number of attempts in structure generation.
            When the number of generated structures is small than ``expect_size``, new rounds of structure generation will be performed.
            The generation will stop until the number of generated structures is more than ``expect_size, `` or the total attempts reach the ``max_attempts``.
            Default ``None``, means ``max_attempts`` equal to parameter ``expect_size``.
        check_distance: bool, optional
            Whether the atomic distance should be checked. default ``True``
        distance_scale_factor : float, optional
            Scale factor to determine the tolerance of atomic distances when distance checking. Unit is Å,
            When ``check_distance`` is ``True``, Any structure has
            all_atomic_distance < (A_atom_covalent_radius + B_atom_covalent_radius) * (1 - distance_scale_factor) will be rejected,
            by default 0.1
        *cfgs: Dict[str, Tuple[str]]
            A tuple with Wyckoff configuration set(s).
            Wyckoff Configuration set is a dict with format like: {"Li": ["a", "c"], "O": ["i"]}.
            Here, the "Li" is an available element symbol and ["a", "c"] is a list
            which contains coresponding Wyckoff letters. For convenience, dict will
            be sorted by keys..
        Yields
        ------
        Tuple[Dict]
            Structure information contains ``spacegroup_mun: int``,
            ``volume: float``, ``lattice: list``, ``wyckoff_letters: list``,
            and ``coords: list``.
        """
        assert expect_size >= 1, 'attempts number must be greater than 1'
        for cfg in cfgs:
            yield cfg, self._cg.gen_many(
                expect_size,
                max_attempts,
                check_distance,
                distance_scale_factor,
                cfg,
            )

    def __repr__(self):
        return f"CrystalGenerator(\
            \n    spacegroup_num={self.spacegroup_num},\
            \n    estimated_volume={self.estimated_volume},\
            \n    estimated_variance={self.estimated_variance},\
            \n    angle_range={self.angle_range},\
            \n    angle_tolerance={self.angle_tolerance},\
            \n    max_attempts_number={self.max_attempts_number},\
            \n    empirical_coords={'...' if self._empirical_coords is not None else None},\
            \n    empirical_coords_variance={self.empirical_coords_variance},\
            \n    empirical_coords_sampling_rate={self.empirical_coords_sampling_rate},\
            \n    empirical_coords_loose_sampling={self.empirical_coords_loose_sampling},\
            \n    verbose={self.verbose}\
            \n    n_jobs={self.n_jobs}\
            \n)"
