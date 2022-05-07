// Copyright 2022 TsumiNa
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

use super::pbc::pbc_all_distances;

use crate::Float;
use crate::RADIUS;
use core::f32::consts::PI;
use ndarray::{arr2, s, Array2, Axis};

pub(crate) fn create_lattice(abc: &Vec<Float>, angles: &Vec<Float>) -> Array2<Float> {
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

pub(crate) fn calculate_volume(abc: &Vec<Float>, angles: &Vec<Float>) -> Float {
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

pub(crate) fn get_multiplied_length(volume: &Float, angles: &Vec<Float>) -> Float {
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

pub(crate) fn get_covalent_radius(element: &Vec<String>) -> Vec<Float> {
    let mut ret: Vec<Float> = Vec::new();
    for e in element.iter() {
        ret.push(RADIUS[e]);
    }
    ret
}

pub(crate) fn check_distance(
    lattice: &Array2<Float>,
    particles: &Array2<Float>,
    elements: &Vec<String>,
    distance_scale_factor: &Float,
) -> bool {
    match pbc_all_distances(lattice, particles) {
        Ok(distance_matrix) => {
            let ii = distance_matrix.shape()[0];
            let radius = get_covalent_radius(elements);
            for i in 0..(ii - 1) {
                for j in (i + 1)..ii {
                    if distance_matrix[[i, j]]
                        < (radius[i] + radius[j]) * (1. - distance_scale_factor)
                    {
                        // if self.verbose {
                        //     println!(
                        //         "[reject] {:>2}[{}] <-> {:>2}[{}] = {:.5} < {:.5} (Å)",
                        //         elements[i],
                        //         wyckoff_letters[i],
                        //         elements[j],
                        //         wyckoff_letters[j],
                        //         distance_matrix[[i, j]],
                        //         (radius[i] + radius[j]) * (1. - distance_scale_factor)
                        //     );
                        // }
                        return false;
                    }
                }
            }

            return true;
        }
        _ => return false,
    }
}

pub(crate) fn lattice_to(lattice: Array2<Float>) -> (Vec<Float>, Vec<Float>) {
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

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ndarray::arr2;

    #[test]
    fn covalent_radius() {
        let elements = vec![
            "Li".to_owned(),
            "P".to_owned(),
            "O".to_owned(),
            "Ti".to_owned(),
            "Pd".to_owned(),
        ];
        assert_eq!(
            get_covalent_radius(&elements),
            vec![1.21, 1.04, 0.64, 1.52, 1.33]
        )
    }

    #[test]
    fn distance_checker() {
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
        assert_eq!(check_distance(&lattice, &particles, &elements, &0.1), false);
    }

    #[test]
    fn lattice_generation() {
        let abc = vec![5.0 as Float, 5.0, 5.0];
        let angles = vec![90.0 as Float, 90.0, 90.0];

        let lattice = create_lattice(&abc, &angles);
        assert_abs_diff_eq!(
            lattice,
            arr2(&[[5.0 as Float, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]]),
            epsilon = 1e-15,
        );
    }

    #[test]
    fn volume() {
        let abc = vec![5.0 as Float, 5.0, 5.0];
        let angles = vec![90.0, 90.0, 90.0];

        let vol = calculate_volume(&abc, &angles);
        assert_abs_diff_eq!(vol, (5.0 as Float).powi(3));
    }

    #[test]
    fn multiplied_length() {
        let vol = (5.0 as Float).powi(3);
        let angles = vec![90.0, 90.0, 90.0];

        assert_abs_diff_eq!(get_multiplied_length(&vol, &angles), (5.0 as Float).powi(3));
    }
}
