# Copyright 2020 TsumiNa. All rights reserved.
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file.

import pytest
from crystallus import WyckoffCfgGenerator as WyRs


def rust_single_thread_no_return(n):
    p = '''(x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
    (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
    (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
    (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
    (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
    (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)'''
    count = WyRs(p, (1 / 3, 1 / 3, 2 / 3), (2 / 3, 2 / 3, 1 / 3)).gen_many(n, True)
    assert count == n


def rust_parallel_no_return(n):
    p = '''(x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
    (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
    (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
    (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
    (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
    (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)'''
    count = WyRs(p, (1 / 3, 1 / 3, 2 / 3), (2 / 3, 2 / 3, 1 / 3)).par_gen_many(n, True)
    assert count == n


def rust_single_thread(n):
    p = '''(x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
    (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
    (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
    (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
    (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
    (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)'''
    count = WyRs(p, (1 / 3, 1 / 3, 2 / 3), (2 / 3, 2 / 3, 1 / 3)).gen_many(n)
    assert len(count) == n


def rust_parallel(n):
    p = '''(x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
    (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
    (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
    (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
    (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
    (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)'''
    count = WyRs(p, (1 / 3, 1 / 3, 2 / 3), (2 / 3, 2 / 3, 1 / 3)).par_gen_many(n)
    assert len(count) == n



sizes = [1_000, 10_000, 100_000, 500_000, 600_000, 700_000, 750_000, 780_000, 790_000, 800_000]
# sizes = [1000, 10_000]
cases = [
    rust_single_thread_no_return, rust_parallel_no_return, rust_single_thread, rust_parallel,
    # py_single_thread_without_numpy, py_single_thread, py_parallel_without_numpy, py_parallel
]


@pytest.mark.parametrize('n', sizes)
@pytest.mark.parametrize('case', cases)
def test_(benchmark, n, case):
    benchmark(case, n)
