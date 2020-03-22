from typing import Dict, List, Tuple
from ._util import dangerwrap
from .crystallus import WyckoffCfgGenerator as _WYG


class WyckoffCfgGenerator(object):

    def __init__(self, *, max_recurrent=1_000, n_jobs=-1, **composition):
        """A generator to probale Wyckoff configurations for the given composition.

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
        return self._wyg.gen_one(spacegroup_num)

    def gen_many(
        self,
        size: int,
        *spacegroup_num: int,
        iterative: bool = False,
    ):
        """Generate probale configurations

        Parameters
        ----------
        size : int
            How many times for trying.
        iterative: bool
            Running like a iterator. Instead of return results until all done,
            Results will be returned when generation of each spacegroup has done.
            By default, False
        spacegroup_num: int
            The spacegroup numbers.

        Returns
        -------
        Dict[int, List[Dict]], List[Dict]
            A collection contains spacegroup number and it's corresponding Wyckoff
            configurations (wy_cfg). If only one spacegroup number was given,
            will only return the list of wy_cfgs, otherwise return in dict with
            spacegroup number as key. wy_cfgs will be formated as
            {<element 1>: (<Wyckoff letter>, <Wyckoff letter>, ...), <element 2>: (...), ...}.
        """
        if iterative:
            for sp_num in spacegroup_num:
                yield sp_num, self._wyg.gen_many(size, sp_num)
        else:
            return self._wyg.gen_many(size, *spacegroup_num)

    def __repr__(self):
        return f"WyckoffCfgGenerator(\
            \n    max_recurrent={self.max_recurrent},\
            \n    n_jobs={self.n_jobs}\
            \n    composition={self.composition}\
            \n)"