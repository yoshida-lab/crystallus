// Copyright 2020 TsumiNa. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

use crate::{Float, SPG_TYPES, WY};
use rand::{thread_rng, Rng};
use std::collections::{BTreeMap, HashMap};
use std::{error, fmt};

#[derive(Debug, Clone)]
pub struct WyckoffCfgGeneratorError(String);

impl fmt::Display for WyckoffCfgGeneratorError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "wyckoff configuration generator error: `{}`", self.0)
    }
}

impl error::Error for WyckoffCfgGeneratorError {}

pub struct WyckoffCfgGenerator<'a> {
    // pool:       [multiplicity, letter, reuse]
    candidate_pool: Vec<(usize, &'a str, &'a bool)>,
    pub max_recurrent: u16,
    scale: u8,
}

impl<'a> WyckoffCfgGenerator<'a> {
    pub fn from_spacegroup_num(
        spacegroup_num: usize,
        max_recurrent: Option<u16>,
    ) -> Result<WyckoffCfgGenerator<'a>, WyckoffCfgGeneratorError> {
        if !(1..=230).contains(&spacegroup_num) {
            return Err(WyckoffCfgGeneratorError(
                "space group number is illegal".to_owned(),
            ));
        }

        // get wyckoff letters and corresponding multiplicity
        //                 [multiplicity, letter, reuse]
        let candidate_pool: Vec<(usize, &'a str, &bool)> = WY[spacegroup_num - 1]
            .iter()
            .map(|(&letter, (multiplicity, reuse, _))| (*multiplicity, letter, reuse))
            .collect();

        let scale = match SPG_TYPES[spacegroup_num - 1] {
            'A' => 2,
            'B' => 2,
            'C' => 2,
            'I' => 2,
            'S' => 3,
            'T' => 3,
            'F' => 4,
            _ => 1,
        };

        Ok(WyckoffCfgGenerator {
            candidate_pool,
            max_recurrent: max_recurrent.unwrap_or(1000),
            scale,
        })
    }

    pub fn gen(
        &self,
        composition: &BTreeMap<&'a str, Float>,
    ) -> Result<BTreeMap<&str, Vec<&str>>, WyckoffCfgGeneratorError> {
        let mut checker = false;
        let mut ret_: BTreeMap<&str, Vec<&str>> = BTreeMap::new(); // such like: {Li: (b, a), O: (d,)}, where 'b', 'a', and 'd' are wyckoff letters
        let mut used: Vec<&str> = Vec::new(); // such like: [b, d]
        let mut composition_: HashMap<&str, Float> = composition
            .iter()
            .map(|(&k, v)| (k, v * self.scale as Float))
            .collect();
        for _ in 0..self.max_recurrent {
            if checker {
                ret_ = BTreeMap::new(); // such like: {Li: (b, a), O: (d,)}, where 'b', 'a', and 'd' are wyckoff letters
                used = Vec::new(); // such like: [b, d]
                composition_ = composition
                    .iter()
                    .map(|(&k, v)| (k, v * self.scale as Float))
                    .collect();
            }

            for (element, num) in composition_.iter_mut() {
                // pool: [multiplicity, letter, reuse]
                let pool: Vec<&(usize, &str, &bool)> = self
                    .candidate_pool
                    .iter()
                    .filter(|(multiplicity, letter, _)| {
                        // only the reuseable and multiplicity smaller than atom number letters
                        !used.contains(letter) && *multiplicity <= *num as usize
                    })
                    .collect();

                let upper = pool.len();
                if upper == 0 {
                    checker = true;
                    continue;
                }
                // randomly select a suitable wyckoff letter
                let n: usize = thread_rng().gen_range(0, pool.len());
                let &(multiplicity, letter, reuse) = pool[n];

                // atom numbers - multiplicity
                *num -= multiplicity as Float;

                // update returns
                ret_.entry(*element).or_insert(vec![]).push(letter);

                // record non-reuseable letters
                if !reuse {
                    used.push(letter);
                }
            }

            composition_.retain(|_, &mut num| num > 0.);
            if composition_.is_empty() {
                return Ok(ret_);
            }
        }

        return Err(WyckoffCfgGeneratorError(
            "reached the max recurrent number".to_owned(),
        ));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::iter::FromIterator;
    #[test]
    fn test_wyckoff_cfg_generator() -> Result<(), WyckoffCfgGeneratorError> {
        let wy = WyckoffCfgGenerator::from_spacegroup_num(30, None)?; // Pnc2
        assert_eq!(wy.scale, 1);
        let wy = WyckoffCfgGenerator::from_spacegroup_num(65, None)?; // Cmmm
        assert_eq!(wy.scale, 2);
        let wy = WyckoffCfgGenerator::from_spacegroup_num(167, None)?; // R-3c
        assert_eq!(wy.scale, 1);
        let wy = WyckoffCfgGenerator::from_spacegroup_num(227, None)?; // Fd-3m
        assert_eq!(wy.scale, 4);
        Ok(())
    }

    #[test]
    fn test_wyckoff_cfg_gen() -> Result<(), WyckoffCfgGeneratorError> {
        let wy = WyckoffCfgGenerator::from_spacegroup_num(167, None)?; // R-3c
        let cfg = wy.gen(&BTreeMap::from_iter(
            vec![("Ca", 2 as Float), ("C", 2.), ("O", 6.)].into_iter(),
        ))?;
        assert!(
            cfg == BTreeMap::from_iter(
                vec![("Ca", vec!["b"]), ("C", vec!["a"]), ("O", vec!["e"])].into_iter(),
            ) || cfg
                == BTreeMap::from_iter(
                    vec![("Ca", vec!["b"]), ("C", vec!["a"]), ("O", vec!["d"])].into_iter(),
                )
                || cfg
                    == BTreeMap::from_iter(
                        vec![("Ca", vec!["a"]), ("C", vec!["b"]), ("O", vec!["e"])].into_iter(),
                    )
                || cfg
                    == BTreeMap::from_iter(
                        vec![("Ca", vec!["a"]), ("C", vec!["b"]), ("O", vec!["d"])].into_iter(),
                    )
        );
        Ok(())
    }
}
