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

mod crystal;
mod options;

use core::f32::consts::PI;
use rand::distributions::{Bernoulli, Distribution};

pub use self::crystal::Crystal;
pub use self::options::CrystalGeneratorOption;
use crate::utils::pbc_all_distances;
use crate::{wyckoff_pos::*, Float, SPG_TYPES, WY};
use error::Error;
use ndarray::{arr2, concatenate, s, Array, Array2, ArrayView2, Axis};
use rand::{thread_rng, Rng};
use rand_distr::{Normal, Uniform};
use std::collections::HashMap;
use std::{error, fmt};

// Covalent radius for element H (Z=1) to Cm (Z=96)
const COVALENT_RADIUS: &'static str = std::include_str!("covalent_radius.json");
lazy_static! {
    static ref RADIUS: HashMap<String, Float> = serde_json::from_str(COVALENT_RADIUS).unwrap();
}

#[derive(Debug, Clone)]
pub struct CrystalGeneratorError(String);

impl fmt::Display for CrystalGeneratorError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "crystal generator error: `{}`", self.0)
    }
}

impl Error for CrystalGeneratorError {}

pub type LatticeFn =
    Box<dyn Fn() -> Result<(Array2<Float>, Float), CrystalGeneratorError> + Send + Sync>;

// #[derive(Debug, Clone, PartialEq)]
pub struct CrystalGenerator {
    pub spacegroup_num: usize,
    pub max_attempts_number: u16,
    pub verbose: bool,
    empirical_coords: HashMap<String, Vec<Vec<Float>>>,
    empirical_coords_variance: Float,
    empirical_coords_sampling_rate: Float,
    empirical_coords_loose_sampling: bool,
    wy_pos_generator: HashMap<String, (usize, WyckoffPos)>,
    lattice_gen: LatticeFn,
}

impl<'a> CrystalGenerator {
    pub fn from_spacegroup_num(
        spacegroup_num: usize,
        estimated_volume: Float,
        estimated_variance: Float,
        options: CrystalGeneratorOption,
    ) -> Result<CrystalGenerator, CrystalGeneratorError> {
        if !(1..=230).contains(&spacegroup_num) {
            return Err(CrystalGeneratorError(
                "space group number is illegal".to_owned(),
            ));
        }
        let verbose = options.verbose;
        let shifts = match &SPG_TYPES[spacegroup_num - 1] {
            'A' => vec![[0., 1. / 2., 1. / 2.]],
            'B' => vec![[1. / 2., 0., 1. / 2.]],
            'C' => vec![[1. / 2., 1. / 2., 0.]],
            'I' => vec![[1. / 2., 1. / 2., 1. / 2.]],
            'S' => vec![[1. / 3., 1. / 3., 2. / 3.], [2. / 3., 2. / 3., 1. / 3.]],
            'T' => vec![[1. / 3., 2. / 3., 1. / 3.], [2. / 3., 1. / 3., 2. / 3.]],
            'F' => vec![
                [0., 1. / 2., 1. / 2.],
                [1. / 2., 0., 1. / 2.],
                [1. / 2., 1. / 2., 0.],
            ],
            _ => vec![],
        };
        let coefficient = (shifts.len() + 1) as Float;
        let estimated_volume = estimated_volume * coefficient;
        let estimated_variance = estimated_variance * coefficient;

        // build wyckoff generators
        let mut wy_pos_generator: HashMap<String, (usize, WyckoffPos)> = HashMap::new();
        WY[spacegroup_num - 1].iter().for_each(|(k, (m, _, p))| {
            wy_pos_generator.insert(
                (*k).to_owned(),
                (
                    *m,
                    WyckoffPos::from_str_and_shifts(p, shifts.clone()).unwrap(),
                ),
            );
        });

        // set value for conditional parameters
        let angle_tolerance = options.angle_tolerance;
        let angle_range = options.angle_range;

        if angle_range.0 * 3. >= 360. {
            return Err(CrystalGeneratorError(
                "angle range is illegal, sum of min bound must smaller than 360 degree".to_owned(),
            ));
        }

        let (low_length, high_length) = (
            estimated_volume.powf(1. / 3.) * 0.5,
            estimated_volume.powf(1. / 3.) * 1.5,
        );
        let length_dist = Uniform::new_inclusive(low_length, high_length);
        let angle_dist = Uniform::new_inclusive(angle_range.0, angle_range.1);
        let volume_dist = Normal::new(estimated_volume, estimated_variance).unwrap();

        let empirical_coords: HashMap<String, Vec<Vec<Float>>> =
            if options.empirical_coords.len() > 0 {
                let mut ret = HashMap::new();
                let coords = options.empirical_coords;
                for (wy, coord) in coords.into_iter() {
                    ret.entry(wy).or_insert(vec![]).push(coord);
                }
                ret
            } else {
                HashMap::new()
            };
        let empirical_coords_variance = options.empirical_coords_variance;
        let empirical_coords_sampling_rate = options.empirical_coords_sampling_rate;
        let empirical_coords_loose_sampling = options.empirical_coords_loose_sampling;

        let max_attempts_number = options.max_attempts_number;

        // generate `lattice generation` function
        // note that closure in Rust is not a available type
        // we have to use Box to wrap and save closure
        let lattice: Array2<Float> =
            Array::from_shape_vec((3, 3), options.lattice).or_else(|e| {
                Err(CrystalGeneratorError(
                    format!("`lattice` is illegal: {}", e).to_owned(),
                ))
            })?;
        let lattice_gen: LatticeFn = if lattice.sum() > 0. {
            // let lattice = options.lattice;
            let (abc, angles) = Self::lattice_to(lattice);
            let vol_ = Self::_volume(&abc, &angles);

            Box::new(move || {
                let vol: Float = thread_rng().sample(volume_dist).abs();
                let abc = abc.clone();
                let angles = angles.clone();
                let vol_ = vol_;

                let ratio = (vol / vol_).powf(1. / 3.);
                let abc = abc.into_iter().map(|x| x * ratio).collect();

                Ok((CrystalGenerator::lattice_from(abc, angles), vol))
            })
        } else {
            match spacegroup_num {
                // Triclinic, α≠β≠γ≠；a≠b≠c
                1..=2 => Box::new(move || {
                    let vol: Float = thread_rng().sample(volume_dist).abs();

                    for _ in 0..max_attempts_number {
                        let angles: Vec<Float> =
                            thread_rng().sample_iter(angle_dist).take(3).collect();
                        // gen angles conditional
                        let sum = angles.iter().sum::<Float>();
                        if angle_tolerance < sum
                            && sum < 360. - angle_tolerance
                            && angles
                                .iter()
                                .filter(|&&x| x > (sum / 2. - angle_tolerance))
                                .count()
                                == 0
                        {
                            let mut abc: Vec<Float> =
                                thread_rng().sample_iter(length_dist).take(2).collect();
                            // c = a*b*c / a*b
                            let c: Float = CrystalGenerator::get_multiplied_length(&vol, &angles)
                                / abc.iter().product::<Float>();

                            if low_length < c && c < high_length {
                                abc.push(c);
                                return Ok((CrystalGenerator::lattice_from(abc, angles), vol));
                            }
                        }
                    }
                    return Err(CrystalGeneratorError(
                        "reached the max attempts (in lattice generation)".to_owned(),
                    ));
                }),
                // Monoclinic, α=γ=90°，β≠90°；a≠b≠c
                3..=15 => Box::new(move || {
                    let vol: Float = thread_rng().sample(volume_dist).abs();

                    // gen angles conditional
                    let mut angles = vec![90 as Float, 0., 90.];
                    for _ in 0..max_attempts_number {
                        let beta = thread_rng().sample(angle_dist);

                        if beta != 90. {
                            angles[1] = beta;
                            let mut abc: Vec<Float> =
                                thread_rng().sample_iter(length_dist).take(2).collect();
                            // c = a*b*c / a*b
                            let c: Float = CrystalGenerator::get_multiplied_length(&vol, &angles)
                                / abc.iter().product::<Float>();

                            if low_length < c && c < high_length {
                                abc.push(c);
                                return Ok((CrystalGenerator::lattice_from(abc, angles), vol));
                            }
                        }
                    }
                    return Err(CrystalGeneratorError(
                        "reached the max attempts (in lattice generation)".to_owned(),
                    ));
                }),
                // Orthorhombic, α=β=γ=90°；a≠b≠c
                16..=74 => Box::new(move || {
                    let vol: Float = thread_rng().sample(volume_dist).abs();

                    // gen angles conditional
                    let angles = vec![90 as Float; 3];
                    for _ in 0..max_attempts_number {
                        let mut abc: Vec<Float> =
                            thread_rng().sample_iter(length_dist).take(2).collect();
                        // c = a*b*c / a*b
                        let c: Float = vol / abc.iter().product::<Float>();

                        if low_length < c && c < high_length {
                            abc.push(c);
                            return Ok((CrystalGenerator::lattice_from(abc, angles), vol));
                        }
                    }
                    return Err(CrystalGeneratorError(
                        "reached the max attempts (in lattice generation)".to_owned(),
                    ));
                }),
                // Tetragonal, α=β=γ=90°；a=b≠c
                75..=142 => Box::new(move || {
                    let vol: Float = thread_rng().sample(volume_dist).abs();

                    // gen angles conditional
                    let angles = vec![90 as Float; 3];
                    for _ in 0..max_attempts_number {
                        let a = thread_rng().sample(length_dist);
                        let c: Float = vol / (a * a);

                        if low_length < c && c < high_length {
                            let abc = vec![a, a, c];
                            return Ok((CrystalGenerator::lattice_from(abc, angles), vol));
                        }
                    }
                    return Err(CrystalGeneratorError(
                        "reached the max attempts (in lattice generation)".to_owned(),
                    ));
                }),
                // Trigonal
                n @ 143..=167 => Box::new(move || {
                    let vol: Float = thread_rng().sample(volume_dist).abs();
                    if ![146, 148, 155, 160, 161, 166, 167].contains(&n) {
                        // α=β=90°，γ=120°；a=b≠c
                        // gen angles conditional
                        let angles = vec![90 as Float, 90., 120.];
                        for _ in 0..max_attempts_number {
                            let a = thread_rng().sample(length_dist);
                            // c = V / a^2 sinγ
                            let sin_gamma = (120 as Float).to_radians().sin();
                            let c: Float = vol / (a * a * sin_gamma);

                            if low_length < c && c < high_length {
                                let abc = vec![a, a, c];
                                return Ok((CrystalGenerator::lattice_from(abc, angles), vol));
                            }
                        }
                        return Err(CrystalGeneratorError(
                            "reached the max attempts (in lattice generation)".to_owned(),
                        ));
                    } else {
                        // α=β=γ<90°；a=b=c
                        for _ in 0..max_attempts_number {
                            // gen angles conditional
                            let angles =
                                vec![thread_rng().sample(Uniform::new(angle_range.0, 90.)); 3];
                            let c: Float = CrystalGenerator::get_multiplied_length(&vol, &angles)
                                .powf(1. / 3.);

                            if low_length < c && c < high_length {
                                let abc = vec![c, c, c];
                                return Ok((CrystalGenerator::lattice_from(abc, angles), vol));
                            }
                        }
                        return Err(CrystalGeneratorError(
                            "reached the max attempts (in lattice generation)".to_owned(),
                        ));
                    }
                }),
                // Hexagonal, α=β=90°，γ=120°；a=b≠c
                168..=194 => Box::new(move || {
                    let vol: Float = thread_rng().sample(volume_dist).abs();
                    // gen angles conditional
                    let angles = vec![90 as Float, 90., 120.];
                    for _ in 0..max_attempts_number {
                        let a = thread_rng().sample(length_dist);
                        // c = V / a^2 sinγ
                        let sin_gamma = (120 as Float).to_radians().sin();
                        let c: Float = vol / (a * a * sin_gamma);

                        if low_length < c && c < high_length {
                            let abc = vec![a, a, c];
                            return Ok((CrystalGenerator::lattice_from(abc, angles), vol));
                        }
                    }
                    return Err(CrystalGeneratorError(
                        "reached the max attempts (in lattice generation)".to_owned(),
                    ));
                }),
                // Cubic, α=β=γ=90°；a=b=c
                195..=230 => Box::new(move || {
                    let vol: Float = thread_rng().sample(volume_dist).abs();

                    // gen angles conditional
                    let angles = vec![90 as Float; 3];
                    let c: Float = vol.powf(1. / 3.);
                    let abc = vec![c, c, c];

                    return Ok((CrystalGenerator::lattice_from(abc, angles), vol));
                }),
                // others
                _ => return Err(CrystalGeneratorError("unknown error".to_owned())),
            }
        };
        Ok(CrystalGenerator {
            spacegroup_num,
            wy_pos_generator,
            empirical_coords,
            empirical_coords_variance,
            empirical_coords_sampling_rate,
            empirical_coords_loose_sampling,
            verbose,
            lattice_gen,
            max_attempts_number,
        })
    }

    #[inline]
    fn get_covalent_radius(element: &Vec<String>) -> Vec<Float> {
        let mut ret: Vec<Float> = Vec::new();
        for e in element.iter() {
            ret.push(RADIUS[e]);
        }
        ret
    }

    #[inline]
    fn lattice_from(abc: Vec<Float>, angles: Vec<Float>) -> Array2<Float> {
        let (a, b, c, alpha, beta, gamma) = (
            abc[0],
            abc[1],
            abc[2],
            angles[0].to_radians(),
            angles[1].to_radians(),
            angles[2].to_radians(),
        );
        let (cos_alpha, cos_beta, cos_gamma) = (alpha.cos(), beta.cos(), gamma.cos());
        let (sin_alpha, sin_beta) = (alpha.sin(), beta.sin());
        let gamma_star = (cos_alpha * cos_beta - cos_gamma) / (sin_alpha * sin_beta);
        // Sometimes rounding errors result in values slightly > 1.
        let gamma_star = gamma_star.min(1.).max(-1.).acos();
        let vector_a = [a * sin_beta, 0.0, a * cos_beta];
        let vector_b = [
            -b * sin_alpha * gamma_star.cos(),
            b * sin_alpha * gamma_star.sin(),
            b * cos_alpha,
        ];
        let vector_c = [0.0, 0.0, c];

        arr2(&[vector_a, vector_b, vector_c])
    }

    #[inline]
    fn lattice_to(lattice: Array2<Float>) -> (Vec<Float>, Vec<Float>) {
        fn abs_cap(val: Float) -> Float {
            val.min(1.).max(-1.)
        }

        let abc = lattice
            .mapv(|a| a.powi(2))
            .sum_axis(Axis(1))
            .mapv(f64::sqrt)
            .into_raw_vec();

        let mut angles = vec![0. as Float, 0., 0.];
        for i in 0..3 {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            let mj = lattice.slice(s![j, ..]);
            let mk = lattice.slice(s![k, ..]);
            angles[i] = abs_cap(mj.dot(&mk) / (abc[j] * abc[k]));
        }
        angles = angles
            .into_iter()
            .map(|x| x.acos() * 180.0 / PI as Float)
            .collect();
        (abc, angles)
    }

    #[inline]
    pub fn gen(
        &self,
        elements: &Vec<String>,
        wyckoff_letters: &Vec<String>,
        check_distance: Option<bool>,
        distance_scale_factor: Option<Float>,
    ) -> Result<Crystal, CrystalGeneratorError> {
        // set default value for options
        let check_distance = check_distance.unwrap_or(true);
        let distance_scale_factor = distance_scale_factor.unwrap_or(0.1);

        // prepare containers
        let mut elements_: Vec<String> = Vec::new();
        let mut wyckoff_letters_: Vec<String> = Vec::new();
        let mut coord_pool: HashMap<String, Vec<(Float, Float, Float)>> = HashMap::new();
        let mut wy_gens_: Vec<(&WyckoffPos, String)> = Vec::new(); //  [(wyckoff_pos_gen, wyckoff_letter)...]

        // init random generator
        let mut rng = thread_rng();
        let mut empirical_coords = self.empirical_coords.clone();
        let bernoulli = Bernoulli::new(self.empirical_coords_sampling_rate)
            .map_err(|e| CrystalGeneratorError(format!("{}", e)))?;

        // iter all wyckoff letters and elements
        for (e, w) in elements.iter().zip(wyckoff_letters) {
            match self.wy_pos_generator.get(w) {
                Some((multiplicity, wy_gen)) => {
                    elements_.append(&mut vec![e.to_string(); *multiplicity]);
                    wyckoff_letters_.append(&mut vec![w.to_string(); *multiplicity]);

                    // determining sampling key
                    // loose: by wyckoff letter -> 'a'
                    // not loose: by element+wyckoff: -> 'Ca:a'
                    let key = if self.empirical_coords_loose_sampling {
                        (*w).clone()
                    } else {
                        [(*e).clone(), (*w).clone()].join(":")
                    };
                    wy_gens_.push((wy_gen, key.clone()));

                    // TODO: can be optimized?
                    if !wy_gen.is_cached() {
                        // sampling from empirical coordinates
                        let coord: (Float, Float, Float) = match empirical_coords.get_mut(&key) {
                            Some(coords) => {
                                // only sampling when coords exists and need sampling
                                if coords.len() > 0 && bernoulli.sample(&mut rng) {
                                    let n: usize = thread_rng().gen_range(0, coords.len());
                                    let tmp = coords.remove(n); // remove used coords template
                                    let (x, y, z) = (tmp[0], tmp[1], tmp[2]);
                                    let dist =
                                        Normal::new(0., self.empirical_coords_variance).unwrap();
                                    // sampling a perturbation
                                    let perturbation: Vec<Float> =
                                        thread_rng().sample_iter(dist).take(3).collect();
                                    let (p1, p2, p3) =
                                        (perturbation[0], perturbation[1], perturbation[2]);
                                    // return coord + perturbation
                                    (x + p1, y + p2, z + p3)
                                } else {
                                    thread_rng().gen()
                                }
                            }
                            None => thread_rng().gen(),
                        };
                        coord_pool.entry(key).or_insert(vec![]).push(coord);
                    }
                }
                None => {
                    return Err(CrystalGeneratorError(format!(
                        "space group {} dose not has wyckoff letter [{}]",
                        self.spacegroup_num, w
                    )))
                }
            }
        }
        // generate lengths and angles, then transfer to a lattice object
        let (lattice, vol) = (self.lattice_gen)()?;
        // let lattice = CrystalGenerator::lattice_from(abc, angles);

        // generate particles for each element, respectively
        let mut all_particles: Vec<Array2<Float>> = Vec::new();
        for (g, key) in wy_gens_.iter() {
            // if wyckoff position is fixed
            // just use the cached values

            if g.is_cached() {
                all_particles.push(g.random_gen());
            } else {
                // try to get empirical coordinate
                let coord = coord_pool
                    .get_mut(key)
                    .ok_or_else(|| CrystalGeneratorError("Unexpected no `key` error".to_owned()))?;
                let n: usize = thread_rng().gen_range(0, coord.len());
                let (x, y, z) = coord.remove(n);
                all_particles.push(g.gen(x, y, z));
            }
        }

        // join all generated particles in their generated order
        let particles = concatenate(
            Axis(0),
            &all_particles
                .iter()
                .map(|p| p.view())
                .collect::<Vec<ArrayView2<Float>>>()[..],
        )
        .unwrap();

        // check distances between all particles,
        // if ok, return generated crystal object
        if !check_distance {
            return Ok(Crystal {
                spacegroup_num: self.spacegroup_num,
                elements: elements_,
                wyckoff_letters: wyckoff_letters_,
                volume: vol,
                lattice,
                particles,
            });
        } else {
            if self.check_distance(
                &lattice,
                &particles,
                &elements_,
                &wyckoff_letters_,
                &distance_scale_factor,
            ) {
                return Ok(Crystal {
                    spacegroup_num: self.spacegroup_num,
                    elements: elements_,
                    wyckoff_letters: wyckoff_letters_,
                    volume: vol,
                    lattice,
                    particles,
                });
            }
            return Err(CrystalGeneratorError(
            "Atomic distance check failed for the randomly generated structure. If you tried many times and still get this error, please try to set `estimated_volume` and/or `distance_scale_factor` bigger. (in crystal structure generation)".to_owned(),
        ));
        }
    }

    #[inline]
    fn _volume(abc: &Vec<Float>, angles: &Vec<Float>) -> Float {
        let (a, b, c) = (abc[0], abc[1], abc[2]);
        let (cos_alpha, cos_beta, cos_gamma) = (
            angles[0].to_radians().cos(),
            angles[1].to_radians().cos(),
            angles[2].to_radians().cos(),
        );
        a * b
            * c
            * (1. + 2. * cos_alpha * cos_beta * cos_gamma
                - cos_alpha.powi(2)
                - cos_beta.powi(2)
                - cos_gamma.powi(2))
            .sqrt()
    }

    #[inline]
    fn get_multiplied_length(volume: &Float, angles: &Vec<Float>) -> Float {
        let (cos_alpha, cos_beta, cos_gamma) = (
            angles[0].to_radians().cos(),
            angles[1].to_radians().cos(),
            angles[2].to_radians().cos(),
        );

        volume
            / (1. + 2. * cos_alpha * cos_beta * cos_gamma
                - cos_alpha.powi(2)
                - cos_beta.powi(2)
                - cos_gamma.powi(2))
            .sqrt()
    }

    fn check_distance(
        &self,
        lattice: &Array2<Float>,
        particles: &Array2<Float>,
        elements: &Vec<String>,
        wyckoff_letters: &Vec<String>,
        distance_scale_factor: &Float,
    ) -> bool {
        match pbc_all_distances(lattice, particles) {
            Ok(distance_matrix) => {
                let ii = distance_matrix.shape()[0];
                let radius = CrystalGenerator::get_covalent_radius(elements);
                for i in 0..(ii - 1) {
                    for j in (i + 1)..ii {
                        if distance_matrix[[i, j]]
                            < (radius[i] + radius[j]) * (1. - distance_scale_factor)
                        {
                            if self.verbose {
                                println!(
                                    "[reject] {:>2}[{}] <-> {:>2}[{}] = {:.5} < {:.5} (Å)",
                                    elements[i],
                                    wyckoff_letters[i],
                                    elements[j],
                                    wyckoff_letters[j],
                                    distance_matrix[[i, j]],
                                    (radius[i] + radius[j]) * (1. - distance_scale_factor)
                                );
                            }
                            return false;
                        }
                    }
                }

                return true;
            }
            _ => return false,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::{assert_abs_diff_eq, assert_ulps_eq};
    use ndarray::arr2;

    #[test]
    fn get_covalent_radius() {
        let elements = vec![
            "Li".to_owned(),
            "P".to_owned(),
            "O".to_owned(),
            "Ti".to_owned(),
            "Pd".to_owned(),
        ];
        assert_eq!(
            CrystalGenerator::get_covalent_radius(&elements),
            vec![1.21, 1.04, 0.64, 1.52, 1.33]
        )
    }

    #[test]
    #[should_panic(
        expected = "called `Result::unwrap()` on an `Err` value: CrystalGeneratorError(\"space group number is illegal\")"
    )]
    fn should_panic_when_using_wrong_spacegroup_number() {
        CrystalGenerator::from_spacegroup_num(300, 100., 10., Default::default()).unwrap();
    }

    #[test]
    fn should_return_space_group_illegal_err() {
        let tmp = CrystalGenerator::from_spacegroup_num(300, 100., 10., Default::default());
        match tmp {
            Err(e) => assert_eq!(
                format!("{}", e),
                "crystal generator error: `space group number is illegal`"
            ),
            _ => assert!(false),
        }
    }

    #[test]
    fn should_return_angle_range_illegal_err() {
        let tmp = CrystalGenerator::from_spacegroup_num(
            1,
            100.,
            10.,
            CrystalGeneratorOption {
                angle_range: (160., 170.),
                ..Default::default()
            },
        );
        match tmp {
            Err(e) => assert_eq!(
                format!("{}", e),
                "crystal generator error: `angle range is illegal, sum of min bound must smaller than 360 degree`"
            ),
            _ => assert!(false),
        }
    }
    #[test]
    #[should_panic(expected = "reached the max attempts (in lattice generation)")]
    fn should_panic_when_reach_max_attempts_number() {
        let tmp = CrystalGenerator::from_spacegroup_num(
            2,
            1000.,
            10.,
            CrystalGeneratorOption {
                angle_range: (119.999, 200.),
                ..Default::default()
            },
        )
        .unwrap();
        (tmp.lattice_gen)().unwrap();
    }

    #[test]
    fn create_crystal_generator_from_spacegroup_with_option() -> Result<(), CrystalGeneratorError> {
        let tmp = CrystalGenerator::from_spacegroup_num(
            2,
            1000.,
            0.,
            CrystalGeneratorOption {
                angle_range: (70., 71.),
                ..Default::default()
            },
        )
        .unwrap();
        let (lattice, vol) = (tmp.lattice_gen)()?;
        let (abc, angles) = CrystalGenerator::lattice_to(lattice);

        assert_ulps_eq!(vol, 1000.);
        assert!((5. < abc[0]) & (abc[0] < 15.));
        assert!((5. < abc[1]) & (abc[1] < 15.));
        assert!((5. < abc[2]) & (abc[2] < 15.));
        assert!((70. < angles[0]) & (angles[0] < 71.));
        assert!((70. < angles[1]) & (angles[1] < 71.));
        assert!((70. < angles[2]) & (angles[2] < 71.));

        Ok(())
    }

    #[test]
    fn crystal_generator_gen_lattice() -> Result<(), CrystalGeneratorError> {
        let tmp =
            CrystalGenerator::from_spacegroup_num(200, 1000., 10., Default::default()).unwrap();
        let (lattice, _vol) = (tmp.lattice_gen)()?;
        let (abc, angles) = CrystalGenerator::lattice_to(lattice);
        assert!(abc.len() == 3);
        assert!(abc.iter().eq(abc.iter()), "{}", true);
        assert!(angles.len() == 3);
        assert!(angles.iter().eq(angles.iter()), "{}", true);
        assert!(9.5 < abc[0] && abc[0] < 10.5);
        assert_abs_diff_eq!(angles[0], 90., epsilon = 1e-6);
        Ok(())
    }

    #[test]
    fn check_distance() {
        let tmp =
            CrystalGenerator::from_spacegroup_num(12, 145.75949096679688, 20., Default::default())
                .unwrap();
        let lattice = arr2(&[
            [14.019043922424316, 0.0, -6.127918936726928e-07],
            [
                -4.818087404601101e-07,
                6.381750106811523,
                -2.7895515586351394e-07,
            ],
            [0.0, 0.0, 9.891742706298828],
        ]);

        let particles = arr2(&[
            [0., 0.69679058, 0.25],
            [0., -0.69679058, 0.75],
            [0.5, 1.19679058, 0.25],
            [0.5, -0.19679058, 0.75],
            [0., 0.13124257, 0.43134677],
            [0., -0.13124257, 0.93134677],
            [0., 0.13124257, 0.06865323],
            [0., -0.13124257, -0.43134677],
            [0.5, 0.63124257, 0.43134677],
            [0.5, 0.36875743, 0.93134677],
            [0.5, 0.63124257, 0.06865323],
            [0.5, 0.36875743, -0.43134677],
            [0., 0.26911515, 0.7328701],
            [0., -0.26911515, 1.2328701],
            [0., 0.26911515, -0.2328701],
            [0., -0.26911515, -0.7328701],
            [0.5, 0.76911515, 0.7328701],
            [0.5, 0.23088485, 1.2328701],
            [0.5, 0.76911515, -0.2328701],
            [0.5, 0.23088485, -0.7328701],
            [0., 0.03196704, 0.25],
            [0., -0.03196704, 0.75],
            [0.5, 0.53196704, 0.25],
            [0.5, 0.46803296, 0.75],
        ]);
        let mut elements: Vec<String> = Vec::new();
        elements.append(&mut vec!["Ti".to_owned(); 8]);
        elements.append(&mut vec!["O".to_owned(); 16]);
        assert_eq!(
            tmp.check_distance(
                &lattice,
                &particles,
                &elements,
                &vec!["a".to_owned(); 24],
                &0.1
            ),
            false
        );
    }

    #[test]
    fn crystal_generate() -> Result<(), CrystalGeneratorError> {
        /*
        space group 63
        ======================================================
        Multiplicity | Wyckoff letter |      Coordinates
        ------------------------------------------------------
             4       |         b      | (0,1/2,0) (0,1/2,1/2)
        ------------------------------------------------------
             4       |         a      | (0,0,0) (0,0,1/2)
        ------------------------------------------------------
        */
        let cg = CrystalGenerator::from_spacegroup_num(63, 1000., 10., Default::default())?;
        let cry = cg.gen(
            &vec!["Li".to_owned(), "P".to_owned()],
            &vec!["a".to_owned(), "b".to_owned()],
            None,
            None,
        )?;
        assert_eq!(
            cry.elements,
            vec![
                "Li".to_owned(),
                "Li".to_owned(),
                "Li".to_owned(),
                "Li".to_owned(),
                "P".to_owned(),
                "P".to_owned(),
                "P".to_owned(),
                "P".to_owned()
            ]
        );
        assert_eq!(cry.particles.shape(), [8, 3]);
        assert_eq!(
            cry.particles,
            arr2(&[
                [0., 0., 0.],    // Li
                [0., 0., 0.5],   // Li
                [0.5, 0.5, 0.],  // Li
                [0.5, 0.5, 0.5], // Li
                [0., 0.5, 0.],   // P
                [0., 0.5, 0.5],  // P
                [0.5, 1., 0.],   // P
                [0.5, 1., 0.5]   // P
            ])
        );
        Ok(())
    }

    #[test]
    fn crystal_generate_with_lattice() -> Result<(), CrystalGeneratorError> {
        let lattice_old = arr2(&[
            [1.53132450e+01, 0.00000000e+00, 9.37665824e-16],
            [-4.66967683e-16, 7.62616100e+00, 4.66967683e-16],
            [0.00000000e+00, 0.00000000e+00, 1.07431550e+01],
        ]);

        let cg = CrystalGenerator::from_spacegroup_num(
            63,
            1000.,
            0.,
            CrystalGeneratorOption {
                lattice: lattice_old.clone().into_raw_vec(),
                ..Default::default()
            },
        )?;
        let cry = cg.gen(
            &vec!["Li".to_owned(), "P".to_owned()],
            &vec!["a".to_owned(), "b".to_owned()],
            None,
            None,
        )?;
        let lattice_new = cry.lattice;
        let (abc_new, angles_new) = CrystalGenerator::lattice_to(lattice_new);
        let (abc_old, angles_old) = CrystalGenerator::lattice_to(lattice_old);

        assert_abs_diff_eq!(cry.volume, 2000., epsilon = 1e-6);
        assert_abs_diff_eq!(
            abc_new[0] / abc_old[0],
            abc_new[1] / abc_old[1],
            epsilon = 1e-6
        );
        assert_abs_diff_eq!(angles_new[0], angles_old[0], epsilon = 1e-5);
        assert_abs_diff_eq!(angles_new[1], angles_old[1], epsilon = 1e-5);
        assert_abs_diff_eq!(angles_new[2], angles_old[2], epsilon = 1e-5);
        assert_abs_diff_eq!(angles_new[0], 90., epsilon = 1e-5);

        Ok(())
    }

    #[test]
    fn crystal_generate_with_no_matched_template() -> Result<(), CrystalGeneratorError> {
        /*
        space group 63
        ======================================================
        Multiplicity | Wyckoff letter |      Coordinates
        ------------------------------------------------------
             4       |         b      | (0,1/2,0) (0,1/2,1/2)
        ------------------------------------------------------
             4       |         a      | (0,0,0) (0,0,1/2)
        ------------------------------------------------------
        */
        let template: Vec<(String, Vec<Float>)> = vec![
            ("c".to_owned(), vec![0.2, 0., 0.0]),
            ("d".to_owned(), vec![0.4, 0., 0.0]),
        ];
        let cg = CrystalGenerator::from_spacegroup_num(
            63,
            1000.,
            10.,
            CrystalGeneratorOption {
                empirical_coords: template,
                empirical_coords_variance: 0.,
                ..Default::default()
            },
        )?;
        let cry = cg.gen(
            &vec!["Li".to_owned(), "P".to_owned()],
            &vec!["a".to_owned(), "b".to_owned()],
            None,
            None,
        )?;
        assert_eq!(
            cry.elements,
            vec![
                "Li".to_owned(),
                "Li".to_owned(),
                "Li".to_owned(),
                "Li".to_owned(),
                "P".to_owned(),
                "P".to_owned(),
                "P".to_owned(),
                "P".to_owned()
            ]
        );
        assert_eq!(cry.particles.shape(), [8, 3]);
        assert_eq!(
            cry.particles,
            arr2(&[
                [0., 0., 0.],    // Li
                [0., 0., 0.5],   // Li
                [0.5, 0.5, 0.],  // Li
                [0.5, 0.5, 0.5], // Li
                [0., 0.5, 0.],   // P
                [0., 0.5, 0.5],  // P
                [0.5, 1., 0.],   // P
                [0.5, 1., 0.5]   // P
            ])
        );
        Ok(())
    }
    #[test]
    fn crystal_generate_with_loose_without_perturbation() -> Result<(), CrystalGeneratorError> {
        /*
        space group 167
        =========================================================================
        Multiplicity | Wyckoff letter |      Coordinates
        -------------------------------------------------------------------------
             6       |         e      | (x,0,1/4) (0,x,1/4) (-x,-x,1/4)
                     |                | (-x,0,3/4) (0,-x,3/4) (x,x,3/4)
        -------------------------------------------------------------------------
             6       |         d      | (1/2,0,0) (0,1/2,0) (1/2,1/2,0)
                     |                | (0,1/2,1/2) (1/2,0,1/2) (1/2,1/2,1/2)
        -------------------------------------------------------------------------
             4       |         c      | (0,0,z) (0,0,-z+1/2) (0,0,-z) (0,0,z+1/2)
        -------------------------------------------------------------------------
             2       |         b      | (0,0,0) (0,0,1/2)
        -------------------------------------------------------------------------
             2       |         a      | (0,0,1/4) (0,0,3/4)
        -------------------------------------------------------------------------
        */
        let template: Vec<(String, Vec<Float>)> = vec![
            ("c".to_owned(), vec![0.2, 0., 0.0]),
            ("e".to_owned(), vec![0.4, 0., 0.0]),
        ];
        let cg = CrystalGenerator::from_spacegroup_num(
            167,
            1000.,
            10.,
            CrystalGeneratorOption {
                empirical_coords: template,
                empirical_coords_variance: 0.,
                ..Default::default()
            },
        )?;
        let cry = cg.gen(
            &vec!["Li".to_owned(), "P".to_owned()],
            &vec!["c".to_owned(), "e".to_owned()],
            Some(false),
            None,
        )?;
        assert_eq!(
            cry.elements,
            vec![
                "Li".to_owned(),
                "Li".to_owned(),
                "Li".to_owned(),
                "Li".to_owned(),
                "P".to_owned(),
                "P".to_owned(),
                "P".to_owned(),
                "P".to_owned(),
                "P".to_owned(),
                "P".to_owned()
            ]
        );
        assert_eq!(cry.particles.shape(), [10, 3]);
        assert_abs_diff_eq!(
            cry.particles,
            arr2(&[
                [0.2, 0.2, 0.2],  // Li
                [0.3, 0.3, 0.3],  // Li
                [0.8, 0.8, 0.8],  // Li
                [0.7, 0.7, 0.7],  // Li
                [0.4, 0.1, 0.25], // P
                [0.25, 0.4, 0.1], // P
                [0.1, 0.25, 0.4], // P
                [0.6, 0.9, 0.75], // P
                [0.75, 0.6, 0.9], // P
                [0.9, 0.75, 0.6]  // P
            ])
        );
        Ok(())
    }
    #[test]
    fn crystal_generate_with_strict_without_perturbation() -> Result<(), CrystalGeneratorError> {
        let template: Vec<(String, Vec<Float>)> = vec![
            ("Li:c".to_owned(), vec![0.2, 0., 0.0]),
            ("P:e".to_owned(), vec![0.4, 0., 0.0]),
        ];
        let cg = CrystalGenerator::from_spacegroup_num(
            167,
            1000.,
            10.,
            CrystalGeneratorOption {
                empirical_coords: template,
                empirical_coords_variance: 0.,
                empirical_coords_loose_sampling: false,
                ..Default::default()
            },
        )?;
        let cry = cg.gen(
            &vec!["Li".to_owned(), "P".to_owned()],
            &vec!["c".to_owned(), "e".to_owned()],
            Some(false),
            None,
        )?;
        assert_eq!(
            cry.elements,
            vec![
                "Li".to_owned(),
                "Li".to_owned(),
                "Li".to_owned(),
                "Li".to_owned(),
                "P".to_owned(),
                "P".to_owned(),
                "P".to_owned(),
                "P".to_owned(),
                "P".to_owned(),
                "P".to_owned()
            ]
        );
        assert_eq!(cry.particles.shape(), [10, 3]);
        assert_abs_diff_eq!(
            cry.particles,
            arr2(&[
                [0.2, 0.2, 0.2],  // Li
                [0.3, 0.3, 0.3],  // Li
                [0.8, 0.8, 0.8],  // Li
                [0.7, 0.7, 0.7],  // Li
                [0.4, 0.1, 0.25], // P
                [0.25, 0.4, 0.1], // P
                [0.1, 0.25, 0.4], // P
                [0.6, 0.9, 0.75], // P
                [0.75, 0.6, 0.9], // P
                [0.9, 0.75, 0.6]  // P
            ])
        );
        Ok(())
    }

    #[test]
    fn crystal_generate_with_template_with_perturbation() -> Result<(), CrystalGeneratorError> {
        let template: Vec<(String, Vec<Float>)> = vec![
            ("c".to_owned(), vec![0.2, 0., 0.0]),
            ("e".to_owned(), vec![0.4, 0., 0.0]),
        ];
        let cg = CrystalGenerator::from_spacegroup_num(
            167,
            1000.,
            10.,
            CrystalGeneratorOption {
                empirical_coords: template,
                empirical_coords_variance: 0.001,
                ..Default::default()
            },
        )?;
        let cry = cg.gen(
            &vec!["Li".to_owned(), "P".to_owned()],
            &vec!["c".to_owned(), "e".to_owned()],
            Some(false),
            None,
        )?;
        assert_eq!(
            cry.elements,
            vec![
                "Li".to_owned(),
                "Li".to_owned(),
                "Li".to_owned(),
                "Li".to_owned(),
                "P".to_owned(),
                "P".to_owned(),
                "P".to_owned(),
                "P".to_owned(),
                "P".to_owned(),
                "P".to_owned()
            ]
        );
        assert_eq!(cry.particles.shape(), [10, 3]);
        assert_abs_diff_eq!(
            cry.particles,
            arr2(&[
                [0.2, 0.2, 0.2],  // Li
                [0.3, 0.3, 0.3],  // Li
                [0.8, 0.8, 0.8],  // Li
                [0.7, 0.7, 0.7],  // Li
                [0.4, 0.1, 0.25], // P
                [0.25, 0.4, 0.1], // P
                [0.1, 0.25, 0.4], // P
                [0.6, 0.9, 0.75], // P
                [0.75, 0.6, 0.9], // P
                [0.9, 0.75, 0.6]  // P
            ]),
            epsilon = 1e-2
        );
        Ok(())
    }
}
