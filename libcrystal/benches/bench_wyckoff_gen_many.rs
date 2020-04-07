#![feature(test)]
extern crate test;
use libcrystal::{Float, WyckoffPos};
use ndarray::Array2;
use rayon::prelude::*;
use test::Bencher;

#[bench]
fn bench_wyckoff_pos_gen_many_single_thread_1000(b: &mut Bencher) {
    let n: i32 = 1000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}

#[bench]
fn bench_wyckoff_pos_gen_many_single_thread_10000(b: &mut Bencher) {
    let n = 10_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}

#[bench]
fn bench_wyckoff_pos_gen_many_single_thread_100000(b: &mut Bencher) {
    let n = 100_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}

#[bench]
fn bench_wyckoff_pos_gen_many_single_thread_500000(b: &mut Bencher) {
    let n = 500_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}
#[bench]
fn bench_wyckoff_pos_gen_many_single_thread_600000(b: &mut Bencher) {
    let n = 600_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}
#[bench]
fn bench_wyckoff_pos_gen_many_single_thread_700000(b: &mut Bencher) {
    let n = 700_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}

#[bench]
fn bench_wyckoff_pos_gen_many_single_thread_750000(b: &mut Bencher) {
    let n = 750_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}
#[bench]
fn bench_wyckoff_pos_gen_many_single_thread_780000(b: &mut Bencher) {
    let n = 780_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}
#[bench]
fn bench_wyckoff_pos_gen_many_single_thread_790000(b: &mut Bencher) {
    let n = 790_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}
#[bench]
fn bench_wyckoff_pos_gen_many_single_thread_800000(b: &mut Bencher) {
    let n = 800_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}

/* parallel
////////////////////////////////////////////////
*/
#[bench]
fn bench_wyckoff_pos_gen_many_parallel_1000(b: &mut Bencher) {
    let n: i32 = 1000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_par_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}

#[bench]
fn bench_wyckoff_pos_gen_many_parallel_10000(b: &mut Bencher) {
    let n = 10_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_par_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}

#[bench]
fn bench_wyckoff_pos_gen_many_parallel_100000(b: &mut Bencher) {
    let n = 100_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_par_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}

#[bench]
fn bench_wyckoff_pos_gen_many_parallel_500000(b: &mut Bencher) {
    let n = 500_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_par_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}
#[bench]
fn bench_wyckoff_pos_gen_many_parallel_600000(b: &mut Bencher) {
    let n = 600_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_par_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}
#[bench]
fn bench_wyckoff_pos_gen_many_parallel_700000(b: &mut Bencher) {
    let n = 700_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_par_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}

#[bench]
fn bench_wyckoff_pos_gen_many_parallel_750000(b: &mut Bencher) {
    let n = 750_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_par_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}
#[bench]
fn bench_wyckoff_pos_gen_many_parallel_780000(b: &mut Bencher) {
    let n = 780_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_par_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}
#[bench]
fn bench_wyckoff_pos_gen_many_parallel_790000(b: &mut Bencher) {
    let n = 790_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_par_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}
#[bench]
fn bench_wyckoff_pos_gen_many_parallel_800000(b: &mut Bencher) {
    let n = 800_000;
    let tmp = WyckoffPos::from_str_and_shifts(
        "
        (x,y,z), (-x,-y,z), (-x,y,-z), (x,-y,-z), (z,x,y), (z,-x,-y),
        (-z,-x,y), (-z,x,-y), (y,z,x), (-y,z,-x), (y,-z,-x), (-y,-z,x),
        (y+1/2,x+1/2,z+1/2), (-y+1/2,-x+1/2,z+1/2), (y+1/2,-x+1/2,-z+1/2),
        (-y+1/2,x+1/2,-z+1/2), (x+1/2,z+1/2,y+1/2), (-x+1/2,z+1/2,-y+1/2),
        (-x+1/2,-z+1/2,y+1/2), (x+1/2,-z+1/2,-y+1/2), (z+1/2,y+1/2,x+1/2),
        (z+1/2,-y+1/2,-x+1/2), (-z+1/2,y+1/2,-x+1/2), (-z+1/2,-y+1/2,x+1/2)
        ",
        vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
    )
    .unwrap();

    b.iter(|| {
        let n = test::black_box(n);
        let _ = (0..n)
            .into_par_iter()
            .map(|_| tmp.random_gen())
            .collect::<Vec<Array2<Float>>>();
    });
}
