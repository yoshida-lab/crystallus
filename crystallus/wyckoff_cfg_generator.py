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

from typing import Dict, List, Tuple
from .crystallus import WyckoffCfgGenerator as _WYG


class WyckoffCfgGenerator(object):

    def __init__(self, *, max_recurrent=1_000, n_jobs=-1, **composition):
        """A generator for possible Wyckoff configuration generation.

        Parameters
        ----------
        max_recurrent : int, optional
            Max recurrent until generate a reasonable structure, by default 5_000
        n_jobs : int, optional
            Number of cpu cores when parallel calculation, by default -1
        composition: Dict
            Composition of compounds in the primitive cell; should be formated
            as {<element symbol>: <ratio in float>}.
        """

        self._wyg = _WYG(max_recurrent=max_recurrent, n_jobs=n_jobs, **composition)
        self._composition = composition

    @property
    def max_recurrent(self):
        return self._wyg.max_recurrent

    @property
    def n_jobs(self):
        return self._wyg.n_jobs

    @property
    def composition(self):
        return self._composition

    @n_jobs.setter
    def n_jobs(self, n):
        self._wyg.n_jobs = n

    def gen_one(self, spacegroup_num: int):
        """Try to generate a possible Wyckoff configuration under the given space group.

        Parameters
        ----------
        spacegroup_num : int
            Space group number.

        Returns
        -------
        Dict
            Wyckoff configuration set, which is a dict with format like:
            {"Li": ["a", "c"], "O": ["i"]}. Here, the "Li" is an available element
            symbol and ["a", "c"] is a list which contains coresponding Wyckoff
            letters. For convenience, dict will be sorted by keys.
        """
        return self._wyg.gen_one(spacegroup_num)

    def gen_many(self, size: int, *spacegroup_num: int):
        """Try to generate possible Wyckoff configuration sets.

        Parameters
        ----------
        size : int
            How many times to try for one space group.
        spacegroup_num: int
            The spacegroup numbers.

        Returns
        -------
        Dict[int, List[Dict]], List[Dict]
            A collection contains spacegroup number and it's corresponding Wyckoff
            configurations (wy_cfg). If only one spacegroup number was given,
            will only return the list of wy_cfgs, otherwise return in dict with
            spacegroup number as key. wy_cfgs will be formated as
            {element 1: [Wyckoff_letter, Wyckoff_letter, ...], element 2: [...], ...}.
        """
        return self._wyg.gen_many(size, *spacegroup_num)

    def gen_many_iter(self, size: int, *spacegroup_num: int):
        """Try to generate possible Wyckoff configuration sets.

        Parameters
        ----------
        size : int
            How many times to try for one space group.
        spacegroup_num: int
            The spacegroup numbers.

        Yields
        ------
        Dict[int, List[Dict]], List[Dict]
            A collection contains spacegroup number and it's corresponding Wyckoff
            configurations (wy_cfg). If only one spacegroup number was given,
            will only return the list of wy_cfgs, otherwise return in dict with
            spacegroup number as key. wy_cfgs will be formated as
            {element 1: [Wyckoff_letter, Wyckoff_letter, ...], element 2: [...], ...}.
        """
        for sp_num in spacegroup_num:
            yield sp_num, self._wyg.gen_many(size, sp_num)

    def __repr__(self):
        return f"WyckoffCfgGenerator(\
            \n    max_recurrent={self.max_recurrent},\
            \n    n_jobs={self.n_jobs}\
            \n    composition={self.composition}\
            \n)"