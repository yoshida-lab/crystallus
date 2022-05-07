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

use super::options::EmpiricalGeneratorOption;

use crate::structure::{Crystal, CrystalGeneratorError, LatticeFn};
use crate::utils::{
    calculate_volume, check_distance as check_distance_, create_lattice, get_multiplied_length,
    lattice_to,
};
use crate::Generator;

use crate::{wyckoff_pos::*, Float, SPG_TYPES, WY};

use ndarray::{concatenate, Array, Array2, ArrayView2, Axis};
use rand::distributions::{Bernoulli, Distribution};
use rand::{thread_rng, Rng};
use rand_distr::{Normal, Uniform};
use std::collections::HashMap;

// #[derive(Debug, Clone, PartialEq)]
pub struct EmpiricalGenerator {
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

impl<'a> Generator<EmpiricalGeneratorOption> for EmpiricalGenerator {
    fn from_spacegroup_num(
        spacegroup_num: usize,
        volume_of_cell: Float,
        options: Option<EmpiricalGeneratorOption>,
    ) -> Result<EmpiricalGenerator, CrystalGeneratorError> {
        if !(1..=230).contains(&spacegroup_num) {
            return Err(CrystalGeneratorError(
                "space group number is illegal".to_owned(),
            ));
        }

        // make options
        let options = match options {
            Some(options) => options,
            None => EmpiricalGeneratorOption::default(),
        };

        let verbose = options.base_option.verbose;
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
        let variance_of_volume = options.base_option.variance_of_volume;
        let multiple = (shifts.len() + 1) as Float;
        let volume_of_cell = volume_of_cell * multiple;
        let variance_of_volume = variance_of_volume * multiple;

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
        let angle_tolerance = options.base_option.angle_tolerance;
        let angle_range = options.base_option.angle_range;

        if angle_range.0 * 3. >= 360. {
            return Err(CrystalGeneratorError(
                "angle range is illegal, sum of min bound must smaller than 360 degree".to_owned(),
            ));
        }

        let (low_length, high_length) = (
            volume_of_cell.powf(1. / 3.) * 0.5,
            volume_of_cell.powf(1. / 3.) * 1.5,
        );
        let length_dist = Uniform::new_inclusive(low_length, high_length);
        let angle_dist = Uniform::new_inclusive(angle_range.0, angle_range.1);
        let volume_dist = Normal::new(volume_of_cell, variance_of_volume).unwrap();

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

        let max_attempts_number = options.base_option.max_attempts_number;

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
            let (abc, angles) = lattice_to(lattice);
            let vol_ = calculate_volume(&abc, &angles);

            Box::new(move || {
                let vol: Float = thread_rng().sample(volume_dist).abs();
                let abc = abc.clone();
                let angles = angles.clone();
                let vol_ = vol_;

                let ratio = (vol / vol_).powf(1. / 3.);
                let abc = abc.into_iter().map(|x| x * ratio).collect();

                Ok((create_lattice(&abc, &angles), vol))
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
                            let c: Float = get_multiplied_length(&vol, &angles)
                                / abc.iter().product::<Float>();

                            if low_length < c && c < high_length {
                                abc.push(c);
                                return Ok((create_lattice(&abc, &angles), vol));
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
                            let c: Float = get_multiplied_length(&vol, &angles)
                                / abc.iter().product::<Float>();

                            if low_length < c && c < high_length {
                                abc.push(c);
                                return Ok((create_lattice(&abc, &angles), vol));
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
                            return Ok((create_lattice(&abc, &angles), vol));
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
                            return Ok((create_lattice(&abc, &angles), vol));
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
                                return Ok((create_lattice(&abc, &angles), vol));
                            }
                        }
                        return Err(CrystalGeneratorError(
                            "reached the max attempts (in lattice generation)".to_owned(),
                        ));
                    } else {
                        // α=β=γ<90°；a=b=c
                        for _ in 0..max_attempts_number {
                            // gen angles conditional
                            let angles = vec![
                                thread_rng().sample(Uniform::new(
                                    angle_range.0,
                                    angle_range.1.min(90.)
                                ));
                                3
                            ];
                            let c: Float = get_multiplied_length(&vol, &angles).powf(1. / 3.);

                            if low_length < c && c < high_length {
                                let abc = vec![c, c, c];
                                return Ok((create_lattice(&abc, &angles), vol));
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
                            return Ok((create_lattice(&abc, &angles), vol));
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

                    return Ok((create_lattice(&abc, &angles), vol));
                }),
                // others
                _ => return Err(CrystalGeneratorError("unknown error".to_owned())),
            }
        };

        // return self
        Ok(Self {
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

    fn gen(
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
                                    let n: usize = thread_rng().gen_range(0..coords.len());
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
        // let lattice = create_lattice(&abc, &angles);

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
                let n: usize = thread_rng().gen_range(0..coord.len());
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
            if check_distance_(&lattice, &particles, &elements_, &distance_scale_factor) {
                return Ok(Crystal {
                    spacegroup_num: self.spacegroup_num,
                    elements: elements_,
                    wyckoff_letters: wyckoff_letters_,
                    volume: vol,
                    lattice,
                    particles,
                });
            }
            return Err(
                CrystalGeneratorError(format!(
                    "Atomic distance check failed for structure: \n\nLattice: {:#.4},\nAtoms: {:#.4}\nVolume: {:.4}\n\nIf you tried many times and still get this error, please try to set `volume_of_cell` and/or `distance_scale_factor` bigger. (in crystal structure generation)",
                    lattice, particles, vol
                ),
            ));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::BaseGeneratorOption;
    use approx::{assert_abs_diff_eq, assert_ulps_eq};
    use ndarray::arr2;

    #[test]
    #[should_panic(
        expected = "called `Result::unwrap()` on an `Err` value: CrystalGeneratorError(\"space group number is illegal\")"
    )]
    fn should_panic_when_using_wrong_spacegroup_number() {
        EmpiricalGenerator::from_spacegroup_num(
            300,
            100.,
            Some(EmpiricalGeneratorOption {
                base_option: BaseGeneratorOption {
                    variance_of_volume: 10.,
                    ..Default::default()
                },
                ..Default::default()
            }),
        )
        .unwrap();
    }

    #[test]
    fn should_return_space_group_illegal_err() {
        let tmp = EmpiricalGenerator::from_spacegroup_num(
            300,
            100.,
            Some(EmpiricalGeneratorOption {
                base_option: BaseGeneratorOption {
                    variance_of_volume: 10.,
                    ..Default::default()
                },
                ..Default::default()
            }),
        );
        match tmp {
            Err(e) => assert_eq!(
                format!("{}", e),
                "CrystalGeneratorError -- `space group number is illegal`"
            ),
            _ => assert!(false),
        }
    }

    #[test]
    fn should_return_angle_range_illegal_err() {
        let tmp = EmpiricalGenerator::from_spacegroup_num(
            1,
            100.,
            Some(EmpiricalGeneratorOption {
                base_option: BaseGeneratorOption {
                    variance_of_volume: 10.,
                    angle_range: (160., 170.),
                    ..Default::default()
                },
                ..Default::default()
            }),
        );
        match tmp {
            Err(e) => assert_eq!(
                format!("{}", e),
                "CrystalGeneratorError -- `angle range is illegal, sum of min bound must smaller than 360 degree`"
            ),
            _ => assert!(false),
        }
    }
    #[test]
    #[should_panic(expected = "reached the max attempts (in lattice generation)")]
    fn should_panic_when_reach_max_attempts_number() {
        let tmp = EmpiricalGenerator::from_spacegroup_num(
            2,
            1000.,
            Some(EmpiricalGeneratorOption {
                base_option: BaseGeneratorOption {
                    variance_of_volume: 10.,
                    angle_range: (119.999, 200.),
                    ..Default::default()
                },
                ..Default::default()
            }),
        )
        .unwrap();
        (tmp.lattice_gen)().unwrap();
    }

    #[test]
    fn create_crystal_generator_from_spacegroup_with_option() -> Result<(), CrystalGeneratorError> {
        let tmp = EmpiricalGenerator::from_spacegroup_num(
            2,
            1000.,
            Some(EmpiricalGeneratorOption {
                base_option: BaseGeneratorOption {
                    variance_of_volume: 0.,
                    angle_range: (70., 71.),
                    ..Default::default()
                },
                ..Default::default()
            }),
        )
        .unwrap();
        let (lattice, vol) = (tmp.lattice_gen)()?;
        let (abc, angles) = lattice_to(lattice);

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
        let tmp = EmpiricalGenerator::from_spacegroup_num(
            200,
            1000.,
            Some(EmpiricalGeneratorOption {
                base_option: BaseGeneratorOption {
                    variance_of_volume: 10.,
                    ..Default::default()
                },
                ..Default::default()
            }),
        )
        .unwrap();
        let (lattice, _vol) = (tmp.lattice_gen)()?;
        let (abc, angles) = lattice_to(lattice);
        assert!(abc.len() == 3);
        assert!(abc.iter().eq(abc.iter()), "{}", true);
        assert!(angles.len() == 3);
        assert!(angles.iter().eq(angles.iter()), "{}", true);
        assert!(9.5 < abc[0] && abc[0] < 10.5);
        assert_abs_diff_eq!(angles[0], 90., epsilon = 1e-5);
        Ok(())
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
        let cg = EmpiricalGenerator::from_spacegroup_num(
            63,
            1000.,
            Some(EmpiricalGeneratorOption {
                base_option: BaseGeneratorOption {
                    variance_of_volume: 10.,
                    ..Default::default()
                },
                ..Default::default()
            }),
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
    fn crystal_generate_with_lattice() -> Result<(), CrystalGeneratorError> {
        let lattice_old = arr2(&[
            [1.53132450e+01, 0.00000000e+00, 9.37665824e-16],
            [-4.66967683e-16, 7.62616100e+00, 4.66967683e-16],
            [0.00000000e+00, 0.00000000e+00, 1.07431550e+01],
        ]);

        let cg = EmpiricalGenerator::from_spacegroup_num(
            63,
            1000.,
            Some(EmpiricalGeneratorOption {
                base_option: BaseGeneratorOption {
                    variance_of_volume: 0.,
                    ..Default::default()
                },
                lattice: lattice_old.clone().into_raw_vec(),
                ..Default::default()
            }),
        )?;
        let cry = cg.gen(
            &vec!["Li".to_owned(), "P".to_owned()],
            &vec!["a".to_owned(), "b".to_owned()],
            None,
            None,
        )?;
        let lattice_new = cry.lattice;
        let (abc_new, angles_new) = lattice_to(lattice_new);
        let (abc_old, angles_old) = lattice_to(lattice_old);

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
        let cg = EmpiricalGenerator::from_spacegroup_num(
            63,
            1000.,
            Some(EmpiricalGeneratorOption {
                base_option: BaseGeneratorOption {
                    variance_of_volume: 10.,
                    ..Default::default()
                },
                empirical_coords: template,
                empirical_coords_variance: 0.,
                ..Default::default()
            }),
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
        let cg = EmpiricalGenerator::from_spacegroup_num(
            167,
            1000.,
            Some(EmpiricalGeneratorOption {
                base_option: BaseGeneratorOption {
                    variance_of_volume: 10.,
                    ..Default::default()
                },
                empirical_coords: template,
                empirical_coords_variance: 0.,
                ..Default::default()
            }),
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
        let cg = EmpiricalGenerator::from_spacegroup_num(
            167,
            1000.,
            Some(EmpiricalGeneratorOption {
                base_option: BaseGeneratorOption {
                    variance_of_volume: 10.,
                    ..Default::default()
                },
                empirical_coords: template,
                empirical_coords_variance: 0.,
                empirical_coords_loose_sampling: false,
                ..Default::default()
            }),
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
        let cg = EmpiricalGenerator::from_spacegroup_num(
            167,
            1000.,
            Some(EmpiricalGeneratorOption {
                base_option: BaseGeneratorOption {
                    variance_of_volume: 10.,
                    ..Default::default()
                },
                empirical_coords: template,
                empirical_coords_variance: 0.001,
                ..Default::default()
            }),
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
