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
