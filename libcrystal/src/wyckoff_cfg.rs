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

pub struct WyckoffCfgGenerator {
    // pool:       [multiplicity, letter, reuse]
    candidate_pool: Vec<(usize, String, bool)>,
    pub max_recurrent: u16,
    scale: u8,
}

impl WyckoffCfgGenerator {
    pub fn from_spacegroup_num(
        spacegroup_num: usize,
        max_recurrent: Option<u16>,
    ) -> Result<WyckoffCfgGenerator, WyckoffCfgGeneratorError> {
        if !(1..=230).contains(&spacegroup_num) {
            return Err(WyckoffCfgGeneratorError(
                "space group number is illegal".to_owned(),
            ));
        }

        // get wyckoff letters and corresponding multiplicity
        //                 [multiplicity, letter, reuse]
        let candidate_pool: Vec<(usize, String, bool)> = WY[spacegroup_num - 1]
            .iter()
            .map(|(letter, (multiplicity, reuse, _))| (*multiplicity, (*letter).clone(), *reuse))
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

    /// Generate Wyckoff configurations for the given composition.
    pub fn gen(
        &self,
        composition: &BTreeMap<String, Float>,
    ) -> Result<BTreeMap<String, Vec<String>>, WyckoffCfgGeneratorError> {
        let mut checker = false;
        let mut ret_: BTreeMap<String, Vec<String>> = BTreeMap::new(); // such like: {Li: (b, a), O: (d,)}, where 'b', 'a', and 'd' are wyckoff letters
        let mut used: Vec<&String> = Vec::new(); // such like: [b, d]
        let mut composition_: HashMap<String, Float> = composition
            .iter()
            .map(|(k, v)| (k.clone(), v * self.scale as Float))
            .collect();
        for _ in 0..self.max_recurrent {
            if checker {
                ret_ = BTreeMap::new(); // such like: {Li: (b, a), O: (d,)}, where 'b', 'a', and 'd' are wyckoff letters
                used = Vec::new(); // such like: [b, d]
                composition_ = composition
                    .iter()
                    .map(|(k, v)| (k.clone(), v * self.scale as Float))
                    .collect();
            }

            for (element, num) in composition_.iter_mut() {
                // pool: [multiplicity, letter, reuse]
                let pool: Vec<&(usize, String, bool)> = self
                    .candidate_pool
                    .iter()
                    .filter(|(multiplicity, letter, _)| {
                        // only the reuseable and multiplicity smaller than atom number letters
                        !used.contains(&letter) && *multiplicity <= *num as usize
                    })
                    .collect();

                let upper = pool.len();
                if upper == 0 {
                    checker = true;
                    continue;
                }
                // randomly select a suitable wyckoff letter
                let n: usize = thread_rng().gen_range(0, upper);
                let (multiplicity, letter, reuse) = &pool[n];

                // record non-reuseable letters
                if !reuse {
                    used.push(letter);
                }
                // atom numbers - multiplicity
                *num -= *multiplicity as Float;

                // update returns
                ret_.entry(element.clone())
                    .or_insert(vec![])
                    .push((*letter).clone());
            }

            composition_.retain(|_, &mut num| num > 0.);
            if composition_.is_empty() {
                // resort wyckoff letters
                // because iterator is lazy, we have to use for statement
                for (_, k) in ret_.iter_mut() {
                    k.sort();
                }
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
            vec![
                ("Ca".to_owned(), 2 as Float),
                ("C".to_owned(), 2.),
                ("O".to_owned(), 6.),
            ]
            .into_iter(),
        ))?;
        assert!(
            cfg == BTreeMap::from_iter(
                vec![
                    ("Ca".to_owned(), vec!["b".to_owned()]),
                    ("C".to_owned(), vec!["a".to_owned()]),
                    ("O".to_owned(), vec!["e".to_owned()])
                ]
                .into_iter(),
            ) || cfg
                == BTreeMap::from_iter(
                    vec![
                        ("Ca".to_owned(), vec!["b".to_owned()]),
                        ("C".to_owned(), vec!["a".to_owned()]),
                        ("O".to_owned(), vec!["d".to_owned()])
                    ]
                    .into_iter(),
                )
                || cfg
                    == BTreeMap::from_iter(
                        vec![
                            ("Ca".to_owned(), vec!["a".to_owned()]),
                            ("C".to_owned(), vec!["b".to_owned()]),
                            ("O".to_owned(), vec!["e".to_owned()])
                        ]
                        .into_iter(),
                    )
                || cfg
                    == BTreeMap::from_iter(
                        vec![
                            ("Ca".to_owned(), vec!["a".to_owned()]),
                            ("C".to_owned(), vec!["b".to_owned()]),
                            ("O".to_owned(), vec!["d".to_owned()])
                        ]
                        .into_iter(),
                    )
        );
        Ok(())
    }
}
