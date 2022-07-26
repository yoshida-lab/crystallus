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

use crate::{Float, SPG_TYPES, WY};
use rand::{distributions::WeightedIndex, seq::SliceRandom, thread_rng, Rng};
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
    pub max_recurrent: u16,
    // pool:       [multiplicity, letter, reuse, probability]
    candidate_pool: Vec<(usize, String, bool, Float)>,
    scale: u8,
}

impl WyckoffCfgGenerator {
    /// Initialize from a space group number.
    ///
    /// # Arguments
    ///
    /// * `spacegroup_num` - Specify the space group number. Should in *1 ~ 230*  
    /// * `max_recurrent` - The Maximum number of retries. Option
    /// * `priority` - Given the sampling priorities for each Wyckoff letter. Option
    ///
    /// ```
    pub fn from_spacegroup_num(
        spacegroup_num: usize,
        max_recurrent: Option<u16>,
        priority: Option<HashMap<String, Float>>,
    ) -> Result<WyckoffCfgGenerator, WyckoffCfgGeneratorError> {
        if !(1..=230).contains(&spacegroup_num) {
            return Err(WyckoffCfgGeneratorError(
                "space group number is illegal".to_owned(),
            ));
        }

        let candidate_pool = match priority {
            // in this case, only set wyckoff letter will be used.
            Some(prior) => {
                // get wyckoff letters and corresponding multiplicity
                //                 [multiplicity, letter, reuse, priority]
                let mut candidate_pool: Vec<(usize, String, bool, Float)> = WY[spacegroup_num - 1]
                    .iter()
                    .map(|(letter, (multiplicity, reuse, _))| {
                        (*multiplicity, (*letter).clone(), *reuse, 0.)
                    })
                    .collect();
                // reset probability using `priority`
                for (_, letter, _, proba) in candidate_pool.iter_mut() {
                    if prior.contains_key(letter) {
                        *proba = prior[letter];
                    }
                }
                candidate_pool
            }
            None => {
                // get wyckoff letters and corresponding multiplicity
                //                 [multiplicity, letter, reuse, probability]
                let candidate_pool: Vec<(usize, String, bool, Float)> = WY[spacegroup_num - 1]
                    .iter()
                    .map(|(letter, (multiplicity, reuse, _))| {
                        (*multiplicity, (*letter).clone(), *reuse, 1.)
                    })
                    .collect();
                candidate_pool
            }
        };

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
    ///
    /// # Arguments
    ///
    /// * `composition` - A composition of compound. For example: {'Li': 2, 'O': 6}
    pub fn gen(
        &self,
        composition: &BTreeMap<String, Float>,
    ) -> Result<BTreeMap<String, Vec<String>>, WyckoffCfgGeneratorError> {
        let mut empty_pool = false; // flag for the case of empty wyckoff pool but still have non-assigned atoms
        let mut ret_: BTreeMap<String, Vec<String>> = BTreeMap::new(); // such like: {Li: (b, a), O: (d,)}, where 'b', 'a', and 'd' are wyckoff letters
        let mut used: Vec<&String> = Vec::new(); // such like: [b, d]
        let mut composition_: Vec<(String, Float)> = composition
            .iter()
            .map(|(k, v)| (k.clone(), v * self.scale as Float))
            .collect();
        for _ in 0..self.max_recurrent {
            if empty_pool {
                ret_ = BTreeMap::new(); // such like: {Li: (b, a), O: (d,)}, where 'b', 'a', and 'd' are wyckoff letters
                used = Vec::new(); // such like: [b, d]
                composition_ = composition
                    .iter()
                    .map(|(k, v)| (k.clone(), v * self.scale as Float))
                    .collect();
            }

            // shuffle elements in composition
            let mut rng = rand::thread_rng();
            composition_.shuffle(&mut rng);

            // get each element and their # atoms in composition
            for (element, num) in composition_.iter_mut() {
                // pool: [multiplicity, letter, reuse, priority]
                let pool: Vec<&(usize, String, bool, Float)> = self
                    .candidate_pool
                    .iter()
                    .filter(|(multiplicity, letter, _, proba)| {
                        // only the reuseable and multiplicity smaller than atom number letters
                        !used.contains(&letter) && *multiplicity <= *num as usize && *proba > 0.
                    })
                    .collect();

                // next try when there are no wyckoff letters for selection
                if pool.len() == 0 {
                    empty_pool = true;
                    continue;
                }
                // random sampling a suitable wyckoff letter
                // form `pool` by their probability
                let dist = WeightedIndex::new(pool.iter().map(|(_, _, _, proba)| *proba)).unwrap();
                let n = thread_rng().sample(dist);
                let (multiplicity, letter, reuse, _) = &pool[n];

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

            composition_.retain(|(_, num)| num > &0.);
            if composition_.is_empty() {
                // resort wyckoff letters
                // because iterator is lazy, we have to use for statement
                for (_, k) in ret_.iter_mut() {
                    k.sort();
                }
                return Ok(ret_);
            }
            let mut rng = rand::thread_rng();
            composition_.shuffle(&mut rng);
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
    fn wyckoff_cfg_generator() -> Result<(), WyckoffCfgGeneratorError> {
        let wy = WyckoffCfgGenerator::from_spacegroup_num(30, None, None)?; // Pnc2
        assert_eq!(wy.scale, 1);
        let wy = WyckoffCfgGenerator::from_spacegroup_num(65, None, None)?; // Cmmm
        assert_eq!(wy.scale, 2);
        let wy = WyckoffCfgGenerator::from_spacegroup_num(167, None, None)?; // R-3c
        assert_eq!(wy.scale, 1);
        let wy = WyckoffCfgGenerator::from_spacegroup_num(227, None, None)?; // Fd-3m
        assert_eq!(wy.scale, 4);
        Ok(())
    }

    #[test]
    fn wyckoff_cfg_gen_without_priority() -> Result<(), WyckoffCfgGeneratorError> {
        let wy = WyckoffCfgGenerator::from_spacegroup_num(167, None, None)?; // R-3c
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

    #[test]
    fn wyckoff_cfg_gen_with_priority() -> Result<(), WyckoffCfgGeneratorError> {
        let prior = HashMap::from_iter(
            vec![
                ("f".to_owned(), 1.),
                ("d".to_owned(), 1.),
                ("c".to_owned(), 1.),
                ("b".to_owned(), 1.),
                ("a".to_owned(), 1.),
                ("e".to_owned(), 0.),
            ]
            .into_iter(),
        );
        let wy = WyckoffCfgGenerator::from_spacegroup_num(167, None, Some(prior))?; // R-3c
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
                    ("O".to_owned(), vec!["d".to_owned()])
                ]
                .into_iter(),
            ) || cfg
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
