// Copyright 2020 TsumiNa. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
use crate::Float;
use ndarray::{arr1, arr2, stack, Array2, ArrayView2, Axis};
use rand::prelude::*;
use regex::Regex;
use std::{error, fmt, str::FromStr};

pub type Coord = [Float; 3];

#[derive(Debug, Clone)]
pub struct ParseWyckoffError<'a>(&'a str, &'a str);

impl fmt::Display for ParseWyckoffError<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "invalid string to be parsed into `{}`, reason: `{}`.",
            self.0, self.1
        )
    }
}

impl error::Error for ParseWyckoffError<'_> {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        // Generic error, underlying cause isn't tracked.
        None
    }
}

#[derive(Debug, PartialEq, Default, Clone)]
struct Coordinate {
    x: Float,
    y: Float,
    z: Float,
    bias: Float,
}

impl Coordinate {
    #[inline]
    fn gen(&self, x: Float, y: Float, z: Float) -> Float {
        self.x * x + self.y * y + self.z * z + self.bias
    }
}

impl FromStr for Coordinate {
    type Err = ParseWyckoffError<'static>;
    #[inline]
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        lazy_static! {
            static ref RE: Regex = Regex::new(
                r"(?P<op>-?\s*(?:\d{1})?)(?P<sym_cons>(?:[xyz])|(?:\d{1}\s*/\s*\d{1})|\d{1})"
            )
            .unwrap();
        }
        let mut ret = Self::default();
        if !RE.is_match(s) {
            return Err(ParseWyckoffError(
                "Coordinate",
                r"the string must be something like `x-2z+1/2`",
            ));
        }
        for caps in RE.captures_iter(s) {
            let op: Float = match &caps["op"] {
                "-" => -1.,
                "" => 1.,
                _ => caps["op"].parse().unwrap(),
            };

            match &caps["sym_cons"] {
                "x" => ret.x = op,
                "y" => ret.y = op,
                "z" => ret.z = op,
                s if s.to_string().contains("/") => {
                    ret.bias = s
                        .matches(char::is_numeric)
                        .map(|s| s.parse().unwrap())
                        .fold(0., |acc, x: Float| if acc == 0. { x } else { acc / x })
                }
                _ => ret.bias = s.parse().unwrap(),
            };
        }
        Ok(ret)
    }
}

#[derive(Debug, PartialEq, Clone)]
struct Particle(Coordinate, Coordinate, Coordinate);

impl Particle {
    #[inline]
    fn gen(&self, x: Float, y: Float, z: Float) -> Coord {
        [
            self.0.gen(x, y, z),
            self.1.gen(x, y, z),
            self.2.gen(x, y, z),
        ]
    }
}

impl FromStr for Particle {
    type Err = ParseWyckoffError<'static>;
    #[inline]
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        lazy_static! {
            static ref RE: Regex = Regex::new(r",\s*").unwrap();
        };
        if !RE.is_match(s) {
            return Err(ParseWyckoffError(
                "Particle",
                r"the string must be something like `x,y,1/2`",
            ));
        };
        let mut ret: Vec<Coordinate> = Vec::new();
        for m in RE.split(s) {
            ret.push(Coordinate::from_str(m)?)
        }
        if ret.len() != 3 {
            return Err(ParseWyckoffError(
                "Particle",
                "particle coordinate must be 3-dim, but got {}",
            ));
        }
        Ok(Self(ret.remove(0), ret.remove(0), ret.remove(0)))
    }
}

#[derive(Debug, PartialEq, Clone)]
pub struct WyckoffPos {
    particles: Vec<Particle>,
    pub _shift: Option<Vec<Coord>>,
    _cached: Option<Array2<Float>>,
}

impl FromStr for WyckoffPos {
    type Err = ParseWyckoffError<'static>;

    #[inline]
    fn from_str(s: &str) -> std::result::Result<Self, <Self as std::str::FromStr>::Err> {
        lazy_static! {
            static ref RE: Regex = Regex::new(r"\((?P<p>.+?)\)").unwrap();
        }
        if !RE.is_match(s) {
            return Err(ParseWyckoffError(
                "WyckoffPos",
                r"the string must be something like `(x,y,1/2)(-y,x-y,1/2)`",
            ));
        };
        let mut ret: Vec<Particle> = Vec::new();
        for p in RE.captures_iter(s) {
            ret.push(Particle::from_str(&p["p"])?);
        }
        Ok(Self {
            particles: ret,
            _cached: None,
            _shift: None,
        })
    }
}

impl WyckoffPos {
    #[inline]
    pub fn from_str_and_shifts(
        s: &str,
        shifts: Vec<Coord>,
    ) -> Result<Self, <Self as std::str::FromStr>::Err> {
        lazy_static! {
            static ref RE: Regex = Regex::new(r"x|y|z").unwrap();
        }
        let mut ret: WyckoffPos = s.parse()?;
        ret._shift = Some(shifts);
        if RE.is_match(s) {
            Ok(ret)
        } else {
            let arr = ret.gen(0., 0., 0.);
            ret._cached = Some(arr);
            Ok(ret)
        }
    }

    #[inline]
    pub fn gen(&self, x: Float, y: Float, z: Float) -> Array2<Float> {
        match &self._cached {
            None => {
                let pos: Vec<Coord> = self.particles.iter().map(|p| p.gen(x, y, z)).collect();
                let pos = arr2(&pos);
                let coords = match &self._shift {
                    None => pos,
                    Some(shifts) => {
                        // build a vec and put raw coords at position 0
                        // auto broadcast arr2 + arr1, and append results into ret
                        let mut ret = vec![pos];
                        ret.append(&mut shifts.iter().map(|m| &ret[0] + &arr1(m)).collect());
                        // stack all rows
                        stack(
                            Axis(0),
                            &ret.iter()
                                .map(|p| p.view())
                                .collect::<Vec<ArrayView2<Float>>>()[..],
                        )
                        .unwrap()
                    }
                };
                // normalization
                coords.mapv_into(|x| match x {
                    x_ if x_ > 1. => x_ - 1.,
                    x_ if x_ < 0. => x_ + 1.,
                    _ => x,
                })
            }
            Some(arr) => return arr.clone(),
        }
    }

    #[inline]
    pub fn random_gen(&self) -> Array2<Float> {
        let (x, y, z) = thread_rng().gen::<(Float, Float, Float)>();
        self.gen(x, y, z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    #[test]
    fn test_coordinate_from_err() {
        assert!(Coordinate::from_str("error patten").is_err());
    }
    #[test]
    fn test_coordinate_from() {
        assert_eq!(
            Coordinate::from_str("-z+y-2x+1/2").unwrap(),
            Coordinate {
                x: -2.,
                y: 1.,
                z: -1.,
                bias: 0.5
            }
        );
        assert_eq!(
            Coordinate::from_str("1").unwrap(),
            Coordinate {
                x: 0.,
                y: 0.,
                z: 0.,
                bias: 1.
            }
        )
    }
    #[test]
    fn test_particle_from_err() {
        assert!(Particle::from_str("error patten").is_err(),)
    }
    #[test]
    fn test_particle_from() {
        assert_eq!(
            Particle::from_str("y-2x+1/2, x-2z+1/2, z-2y+1/2").unwrap(),
            Particle(
                Coordinate::from_str("y-2x+1/2").unwrap(),
                Coordinate::from_str("x-2z+1/2").unwrap(),
                Coordinate::from_str("z-2y+1/2").unwrap(),
            )
        )
    }
    #[test]
    fn test_wyckoff_pos_from_err() {
        assert!(WyckoffPos::from_str("error patten").is_err())
    }
    #[test]
    fn test_wyckoff_pos_from() {
        assert_eq!(
            WyckoffPos::from_str("(x,y,1/2)(-y,x-y,1/2)(-x+y,-x,1/2)").unwrap(),
            WyckoffPos {
                particles: vec![
                    Particle::from_str("x,y,1/2").unwrap(),
                    Particle::from_str("-y,x-y,1/2").unwrap(),
                    Particle::from_str("-x+y,-x,1/2").unwrap(),
                ],
                _shift: None,
                _cached: None,
            }
        )
    }
    #[test]
    fn test_wyckoff_pos_from_and_shift() {
        assert_eq!(
            WyckoffPos::from_str_and_shifts("(x,y,1/2); (-y,x-y,1/2)", vec![[0.1, 0.1, 0.1,]])
                .unwrap(),
            WyckoffPos {
                particles: vec![
                    Particle::from_str("x,y,1/2").unwrap(),
                    Particle::from_str("-y,x-y,1/2").unwrap(),
                ],
                _shift: Some(vec![[0.1, 0.1, 0.1,]]),
                _cached: None,
            }
        )
    }
    #[test]
    fn test_wyckoff_pos_can_be_cached() {
        assert_eq!(
            WyckoffPos::from_str_and_shifts("(1,1,1/2); (3/4,0,1/2)", vec![])
                .unwrap()
                ._cached,
            Some(arr2(&[[1., 1., 0.5], [0.75, 0., 0.5]])),
        )
    }
    #[test]
    fn test_wyckoff_pos_gen() {
        assert_eq!(
            WyckoffPos::from_str("(x,y,1/2),(-y,x-y,1/2),(-x+y,-x,1/2),(-x,-y,1/2)")
                .unwrap()
                .gen(1., 1., 1.),
            arr2(&[[1., 1., 0.5], [0., 0., 0.5], [0., 0., 0.5], [0., 0., 0.5],])
        );
    }
    #[test]
    fn test_wyckoff_pos_random_gen() {
        let wy = WyckoffPos::from_str("(x,y,1/2),(-y,x-y,1/2)").unwrap();
        assert_ne!(wy.random_gen(), wy.random_gen())
    }
    #[test]
    fn test_wyckoff_pos_gen_with_cached() {
        assert_eq!(
            WyckoffPos::from_str("(1,1,1/2); (3/4,0,1/2)")
                .unwrap()
                .gen(0.1, 0.2, 0.3),
            arr2(&[[1., 1., 0.5], [0.75, 0., 0.5]])
        )
    }
    #[test]
    fn test_gen_with_shift() {
        let tmp = WyckoffPos::from_str_and_shifts(
            "(x,y,1/2),	(-y,x-y,1/2);	(-x+y,-x,1/2)	(-x,-y,1/2)",
            vec![[0.3, 0.3, 0.3]],
        )
        .unwrap()
        .gen(1., 1., 1.);
        assert_abs_diff_eq!(
            tmp,
            arr2(&[
                [1., 1., 0.5],
                [0., 0., 0.5],
                [0., 0., 0.5],
                [0., 0., 0.5],
                [0.3, 0.3, 0.5 + 0.3],
                [0.3, 0.3, 0.5 + 0.3],
                [0.3, 0.3, 0.5 + 0.3],
                [0.3, 0.3, 0.5 + 0.3],
            ])
        )
    }
    #[test]
    fn test_ndarray_to_vec() {
        let a = arr2(&[[1, 2], [3, 4]]);
        assert_eq!(a.into_raw_vec(), [1, 2, 3, 4])
    }
}
