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

use crate::Float;
use ndarray::{arr2, Array, Array2};
use ndarray_linalg::{error::LinalgError, InverseInto};
use std::cmp;

#[inline]
fn _dot_2d_mod(a: &Array2<Float>, b: &Array2<Float>) -> Array2<Float> {
    let ii = a.shape()[0]; // 27
    let jj = b.shape()[1]; // 3
    let kk = a.shape()[1]; // 3

    let mut o = Array::zeros(a.raw_dim());
    for j in 0..jj {
        for k in 0..kk {
            for i in 0..ii {
                let mut tmp = a[[i, k]] % 1.;
                if tmp < 0. {
                    tmp += 1.;
                }
                o[[i, j]] += tmp * b[[k, j]];
            }
        }
    }

    o
}

/// Returns the shortest vectors between two lists of coordinates taking into
/// account periodic boundary conditions and the lattice.
///
/// Args:
///    lattice: lattice to use
///    frac_coords: Fractional coordinates. e.g. [[1.1, 1.2, 4.3], [0.5, 0.6, 0.7]].
///
/// Returns:
///    Distances between each site.
#[inline]
pub(crate) fn pbc_all_distances(
    lattice: &Array2<Float>,
    frac_coords: &Array2<Float>,
) -> Result<Array2<Float>, LinalgError> {
    lazy_static! {
        static ref IMAGES: Array2<Float> = arr2(&[
            [-1., -1., -1.],
            [-1., -1., 0.],
            [-1., -1., 1.],
            [-1., 0., -1.],
            [-1., 0., 0.],
            [-1., 0., 1.],
            [-1., 1., -1.],
            [-1., 1., 0.],
            [-1., 1., 1.],
            [0., -1., -1.],
            [0., -1., 0.],
            [0., -1., 1.],
            [0., 0., -1.],
            [0., 0., 0.],
            [0., 0., 1.],
            [0., 1., -1.],
            [0., 1., 0.],
            [0., 1., 1.],
            [1., -1., -1.],
            [1., -1., 0.],
            [1., -1., 1.],
            [1., 0., -1.],
            [1., 0., 0.],
            [1., 0., 1.],
            [1., 1., -1.],
            [1., 1., 0.],
            [1., 1., 1.]
        ]);
    }

    // create images, 2d array of all length 3 combinations of [-1,0,1]
    let (lll_matrix, mapping) = lll_reduce(lattice, None);
    let lll_inv = mapping.inv_into()?;

    let ii = frac_coords.shape()[0];
    // get inverse coords
    let fc = frac_coords.dot(&lll_inv);

    let cart_f = _dot_2d_mod(&fc, &lll_matrix);
    // let cart_im = IMAGES.dot(lattice);
    let cart_im = IMAGES.dot(&lll_matrix);

    let mut distances = Array2::zeros((ii, ii));
    let mut pre_im = [0 as Float; 3];
    let mut best: Float;
    let mut da: Float;
    let mut db: Float;
    let mut dc: Float;
    let mut d: Float;

    for i in 0..(ii - 1) {
        for j in (i + 1)..ii {
            for l in 0..3 {
                pre_im[l] = cart_f[[j, l]] - cart_f[[i, l]];
            }
            best = 1e100;
            for k in 0..27 {
                // compilers have a hard time unrolling this
                da = pre_im[0] + cart_im[[k, 0]];
                db = pre_im[1] + cart_im[[k, 1]];
                dc = pre_im[2] + cart_im[[k, 2]];
                d = da * da + db * db + dc * dc;
                if d < best {
                    best = d;
                }
            }

            distances[[i, j]] = best.sqrt();
        }
    }

    Ok(distances)
}

#[inline]
fn _gram_schmidt_process(basis: &Array2<Float>) -> Array2<Float> {
    let mut basis = basis.clone();
    for i in 1..basis.shape()[0] {
        for j in 0..i {
            let b_i = basis.row(i);
            let b_j = basis.row(j);
            let alpha = b_j.dot(&b_i) / b_j.dot(&b_j);
            let tmp = &b_i - &(&b_j * alpha).view();
            basis.row_mut(i).assign(&tmp);
        }
    }
    basis
}

/// Lattice reduction using the original Lenstra-Lenstra-Lovasz algorithm
///
/// This function uses platform floating-point numbers (IEEE 754) for all
/// arithmetic operations. The value of `delta` is set to 0.75.
///
///   - `basis`: A generating matrix for the lattice
///
#[inline]
pub(crate) fn lll_reduce(
    basis: &Array2<Float>,
    delta: Option<Float>,
) -> (Array2<Float>, Array2<Float>) {
    // Parameter delta in the Lovasz condition
    let delta = delta.unwrap_or(0.75);
    let mut basis = basis.clone();
    let mut ortho = _gram_schmidt_process(&basis);
    let n = basis.shape()[0];
    let mut k = 1;
    let mut mapping: Array2<Float> = Array::eye(3);

    while k < n {
        for j in (0..k).rev() {
            let o_j = ortho.row(j);
            let b_k = basis.row(k);
            // A proj_coff B := A.dot(B) / A.dot(A)
            let mu_kj = o_j.dot(&b_k) / o_j.dot(&o_j);
            if mu_kj.abs() > 0.5 {
                let round_mu_kj = mu_kj.round();

                let b_k = basis.row(k);
                let b_j = basis.row(j);
                let tmp = &b_k - &(&b_j * round_mu_kj).view();
                basis.row_mut(k).assign(&tmp);

                let m_k = mapping.row(k);
                let m_j = mapping.row(j);
                let tmp = &m_k - &(&m_j * round_mu_kj).view();
                mapping.row_mut(k).assign(&tmp);

                ortho = _gram_schmidt_process(&basis);
            }
        }

        let o_k = ortho.row(k);
        let sdot_o_k = o_k.dot(&o_k);
        let o_k_1 = ortho.row(k - 1);
        let b_k = basis.row(k);
        let sdot_o_k_1 = o_k_1.dot(&o_k_1);
        // A proj_coff B := A.dot(B) / A.dot(A)
        let mu_kk_1 = &o_k_1.dot(&b_k) / sdot_o_k_1;

        if sdot_o_k >= sdot_o_k_1 * (delta - mu_kk_1 * mu_kk_1) {
            k += 1;
        } else {
            for i in 0..n {
                basis.swap([k, i], [k - 1, i]);
                mapping.swap([k, i], [k - 1, i]);
            }
            ortho = _gram_schmidt_process(&basis);
            k = cmp::max(k - 1, 1);
        }
    }
    (basis, mapping)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use ndarray::arr2;
    use ndarray_linalg::{error::LinalgError, InverseInto};
    #[test]
    fn gram_schmidt_process() {
        let ortho =
            _gram_schmidt_process(&arr2(&[[1 as Float, 1., 1.], [-1., 0., 2.], [3., 5., 6.]]));
        // [[1, 1, 1], [-4/3, -1/3, 5/3], [-3/7, 9/14, -3/14]]
        assert_abs_diff_eq!(
            ortho,
            arr2(&[
                [1 as Float, 1., 1.],
                [-4. / 3., -1. / 3., 5. / 3.],
                [-3. / 7., 9. / 14., -3. / 14.]
            ]),
            epsilon = 1e-7
        )
    }
    #[test]
    fn lll_reduce_() -> Result<(), LinalgError> {
        // lattice and frac coords are generated by this generator under space group 63,
        // and the following values are validated by pymatgen
        let lattice = arr2(&[
            [14.019043922424316, 0.0, -6.127918936726928e-07],
            [
                -4.818087404601101e-07,
                6.381750106811523,
                -2.7895515586351394e-07,
            ],
            [0.0, 0.0, 9.891742706298828],
        ]);
        let (basis, mapping) = lll_reduce(&lattice, None);
        assert_abs_diff_eq!(
            basis,
            arr2(&[
                [-4.81808740e-07, 6.38175011e+00, -2.78955156e-07],
                [0.00000000e+00, 0.00000000e+00, 9.89174271e+00],
                [1.40190439e+01, 0.00000000e+00, -6.12791894e-07]
            ]),
            epsilon = 1e-7
        );
        assert_eq!(mapping, arr2(&[[0., 1., 0.], [0., 0., 1.], [1., 0., 0.]]));

        let lll_inv = mapping.inv_into()?;
        let frac_coords = arr2(&[
            [0.0, 0.6967905759811401, 0.25],
            [0.0, -0.6967905759811401, 0.75],
            [0.5, 1.1967905759811401, 0.25],
        ]);
        assert_abs_diff_eq!(
            frac_coords.dot(&lll_inv),
            arr2(&[
                [0.69679058, 0.25, 0.],
                [-0.69679058, 0.75, 0.],
                [1.19679058, 0.25, 0.5]
            ]),
            epsilon = 1e-7
        );
        Ok(())
    }
    #[test]
    fn dot_mod() -> Result<(), LinalgError> {
        // lattice and frac coords are generated by this generator under space group 63,
        // and the following values are validated by pymatgen
        let lattice = arr2(&[
            [14.019043922424316, 0.0, -6.127918936726928e-07],
            [
                -4.818087404601101e-07,
                6.381750106811523,
                -2.7895515586351394e-07,
            ],
            [0.0, 0.0, 9.891742706298828],
        ]);

        let frac_coords = arr2(&[
            [0.0, 0.6967905759811401, 0.25],
            [0.0, -0.6967905759811401, 0.75],
            [0.5, 1.1967905759811401, 0.25],
        ]);
        assert_abs_diff_eq!(
            _dot_2d_mod(&lattice, &frac_coords),
            arr2(&[
                [0.49999969, 1.21005947, 0.25476083],
                [0.49999986, 1.62758061, 0.78631239],
                [0.44587135, 1.06722927, 0.22293568]
            ]),
            epsilon = 1e-5
        );
        Ok(())
    }
    #[test]
    fn pbc_all_distances_() -> Result<(), LinalgError> {
        // lattice and frac coords are generated by this generator under space group 63,
        // and the following values are validated by pymatgen
        let lattice = arr2(&[
            [14.019043922424316, 0.0, -6.127918936726928e-07],
            [
                -4.818087404601101e-07,
                6.381750106811523,
                -2.7895515586351394e-07,
            ],
            [0.0, 0.0, 9.891742706298828],
        ]);

        let frac_coords = arr2(&[
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
        let ret: Vec<Float> = vec![
            0., 5.54711305, 7.70162829, 8.60559507, 3.30226544, 3.33762267, 3.30226557, 3.33762258,
            7.24749803, 7.96562324, 7.24749796, 7.96562355, 5.5012207, 0.27577709, 5.5012205,
            0.27577708, 8.49474441, 7.6159427, 8.49474449, 7.61594272, 2.13901241, 5.24003888,
            7.0880048, 8.70208414, 0., 0., 8.60559507, 7.70162829, 3.33762258, 3.30226557,
            3.33762267, 3.30226544, 7.96562355, 7.24749796, 7.96562324, 7.24749803, 0.27577708,
            5.5012205, 0.27577709, 5.5012207, 7.61594272, 8.49474449, 7.6159427, 8.49474441,
            5.24003888, 2.13901241, 8.70208414, 7.0880048, 0., 0., 0., 5.54711305, 7.24749803,
            7.96562324, 7.24749796, 7.96562355, 3.30226544, 3.33762267, 3.30226557, 3.33762258,
            8.49474441, 7.6159427, 8.49474449, 7.61594272, 5.5012207, 0.27577709, 5.5012205,
            0.27577708, 7.0880048, 8.70208414, 2.13901241, 5.24003888, 0., 0., 0., 0., 7.96562355,
            7.24749796, 7.96562324, 7.24749803, 3.33762258, 3.30226557, 3.33762267, 3.30226544,
            7.61594272, 8.49474449, 7.6159427, 8.49474441, 0.27577708, 5.5012205, 0.27577709,
            5.5012207, 8.70208414, 7.0880048, 5.24003888, 2.13901241, 0., 0., 0., 0., 0.,
            5.22184369, 3.58767118, 2.15655205, 7.70162829, 8.71163395, 8.49626151, 7.2990159,
            3.10966531, 3.22217415, 3.43604418, 3.02763663, 7.96052764, 7.30699839, 8.09360129,
            7.22332354, 1.90242836, 3.31966611, 7.6740567, 7.9804928, 0., 0., 0., 0., 0., 0.,
            2.15655196, 3.58767118, 8.71163395, 7.70162829, 7.29901604, 8.49626151, 3.22217429,
            3.10966539, 3.02763675, 3.43604426, 7.30699828, 7.96052734, 7.22332343, 8.09360096,
            3.31966603, 1.90242842, 7.98049312, 7.6740565, 0., 0., 0., 0., 0., 0., 0., 5.22184369,
            8.49626151, 7.29901604, 7.70162829, 8.71163395, 3.43604426, 3.02763675, 3.10966539,
            3.22217429, 8.09360096, 7.22332343, 7.96052734, 7.30699828, 1.90242842, 3.31966603,
            7.6740565, 7.98049312, 0., 0., 0., 0., 0., 0., 0., 0., 7.2990159, 8.49626151,
            8.71163395, 7.70162829, 3.02763663, 3.43604418, 3.22217415, 3.10966531, 7.22332354,
            8.09360129, 7.30699839, 7.96052764, 3.31966611, 1.90242836, 7.9804928, 7.6740567, 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 5.22184369, 3.58767118, 2.15655205, 7.96052764,
            7.30699839, 8.09360129, 7.22332354, 3.10966531, 3.22217415, 3.43604418, 3.02763663,
            7.6740567, 7.9804928, 1.90242836, 3.31966611, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            2.15655196, 3.58767118, 7.30699828, 7.96052734, 7.22332343, 8.09360096, 3.22217429,
            3.10966539, 3.02763675, 3.43604426, 7.98049312, 7.6740565, 3.31966603, 1.90242842, 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 5.22184369, 8.09360096, 7.22332343, 7.96052734,
            7.30699828, 3.43604426, 3.02763675, 3.10966539, 3.22217429, 7.6740565, 7.98049312,
            1.90242842, 3.31966603, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 7.22332354,
            8.09360129, 7.30699839, 7.96052764, 3.02763663, 3.43604418, 3.22217415, 3.10966531,
            7.9804928, 7.6740567, 3.31966611, 1.90242836, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 5.75724369, 0.33888913, 5.46886635, 7.70162829, 8.58222364, 7.70908061,
            8.39149597, 5.01045829, 1.92888821, 8.64647362, 7.12555932, 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 5.46886613, 0.33888913, 8.58222364, 7.70162829, 8.39149601,
            7.70908061, 1.92888819, 5.01045842, 7.12555933, 8.64647334, 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 5.75724369, 7.70908061, 8.39149601, 7.70162829,
            8.58222364, 5.01045842, 1.92888819, 8.64647334, 7.12555933, 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 8.39149597, 7.70908061, 8.58222364, 7.70162829,
            1.92888821, 5.01045829, 7.12555932, 8.64647362, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 5.75724369, 0.33888913, 5.46886635, 8.64647362, 7.12555932,
            5.01045829, 1.92888821, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 5.46886613, 0.33888913, 7.12555933, 8.64647334, 1.92888819, 5.01045842, 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 5.75724369,
            8.64647334, 7.12555933, 5.01045842, 1.92888819, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 7.12555932, 8.64647362, 1.92888821, 5.01045829,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            4.96267231, 7.70162829, 9.0188339, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 9.0188339, 7.70162829, 0., 0., 0., 0., 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 4.96267231, 0., 0., 0.,
            0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
        ];
        assert_abs_diff_eq!(
            pbc_all_distances(&lattice, &frac_coords)?,
            Array::from_shape_vec((24, 24), ret).unwrap(),
            epsilon = 1e-7
        );
        Ok(())
    }
}
