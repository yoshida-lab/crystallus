// Copyright 2021 TsumiNa
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#![feature(test)]
extern crate test;
use libcrystal::WyckoffPos;
use test::Bencher;

#[bench]
fn bench_wyckoff_pos_gen_one_without_shift(b: &mut Bencher) {
    let tmp: WyckoffPos = "(x,y,1/2)	(-y,x-y,1/2)	(-x+y,-x,1/2)	(-x,-y,1/2)"
        .parse()
        .unwrap();
    b.iter(|| tmp.gen(1., 1., 1.));
}

#[bench]
fn bench_wyckoff_pos_gen_one_with_cached(b: &mut Bencher) {
    let tmp: WyckoffPos = "(1,1,1/2); (3/4,0,1/2)".parse().unwrap();
    b.iter(|| tmp.gen(1., 1., 1.));
}
#[bench]
fn bench_wyckoff_pos_gen_one_with_cached_and_shift(b: &mut Bencher) {
    let tmp: WyckoffPos =
        WyckoffPos::from_str_and_shifts("(1,1,1/2); (3/4,0,1/2)", vec![[0.2, 0.2, 0.2]]).unwrap();
    b.iter(|| tmp.gen(1., 1., 1.));
}
#[bench]
fn bench_wyckoff_pos_gen_one_without_cached(b: &mut Bencher) {
    let tmp: WyckoffPos = "(1,1,x); (3/4,0,y)".parse().unwrap();
    b.iter(|| tmp.gen(1., 1., 1.));
}
#[bench]
fn bench_wyckoff_pos_gen_with_shift(b: &mut Bencher) {
    let tmp = WyckoffPos::from_str_and_shifts(
        "
    (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y), (-z,-x,y), (-z,x,-y), 
    (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x), (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), 
    (y+1/2,-x+1/2,-z+1/2), (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2), 
    (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2), (z+1/2,-y+1/2,-x+1/2), 
    (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)",
        vec![[0.2, 0.2, 0.2]],
    )
    .unwrap();
    b.iter(|| tmp.gen(0.3, 0.4, 0.5));
}
