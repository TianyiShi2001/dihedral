//! # Dihedral
//!
//! [![crates.io](https://img.shields.io/crates/d/dihedral.svg)](https://crates.io/crates/dihedral)
//! [![crates.io](https://img.shields.io/crates/v/dihedral.svg)](https://crates.io/crates/dihedral)
//! [![crates.io](https://img.shields.io/crates/l/dihedral.svg)](https://crates.io/crates/dihedral)
//!
//!  This crate provides functions for working with dihedral angles. Currently, there are two functions:
//!
//! - `dihedral` calculates the dihedral angle in the range -π to π in accordance with biochemistry textbooks (see also: https://en.wikipedia.org/wiki/Dihedral_angle#In_stereochemistry)
//! - `dihedral_unsigned` ignores the direction of rotation and outputs the angle within the range 0 to π. This function is faster than the above signed version.
//!
//! If you want to use `f32` instead of `f64` for calculation, you can add `dihedral = {version = "*", features = ["f32"]}` to your `Cargo.toml`.
//!
//! # References
//!
//! - https://math.stackexchange.com/a/47084
//! - https://en.wikipedia.org/wiki/Dihedral_angle#In_stereochemistry

#[cfg(not(feature = "f32"))]
type Float = f64;
#[cfg(feature = "f32")]
type Float = f32;
type Vector3D = [Float; 3];
type Point3D = [Float; 3];

/// Calculates the dihedral angle, in the range -π to π, of the four ordered coordinates.
///
/// This is in accordance with the definition of dihedral angle in biochemistry textbooks.
///
/// # Examples
///
/// ```
/// use dihedral::dihedral;
///
/// let P0 = [24.969, 13.428, 30.692]; // N
/// let P1 = [24.044, 12.661, 29.808]; // CA
/// let P2 = [22.785, 13.482, 29.543]; // C
/// let P3 = [21.951, 13.670, 30.431]; // O
/// let P4 = [23.672, 11.328, 30.466]; // CB
/// let P5 = [22.881, 10.326, 29.620]; // CG
/// let P6 = [23.691, 9.935, 28.389]; // CD1
/// let P7 = [22.557, 9.096, 30.459]; // CD2
///
/// assert!((dihedral([P0, P1, P2, P3]).to_degrees() - (-71.21515)).abs() < 1E-2);
/// assert!((dihedral([P0, P1, P4, P5]).to_degrees() - (-171.94319)).abs() < 1E-2);
/// assert!((dihedral([P1, P4, P5, P6]).to_degrees() - (60.82226)).abs() < 1E-2);
/// assert!((dihedral([P1, P4, P5, P7]).to_degrees() - (-177.63641)).abs() < 1E-2);
/// ```
#[inline]
pub fn dihedral([a, b, c, d]: [Point3D; 4]) -> Float {
    let (a, b, c) = (v(a, b), v(b, c), v(c, d));
    let (r, s) = (u(cross(a, b)), u(cross(b, c)));
    let t = cross(r, u(b));
    let (x, y) = (dot(r, s), dot(s, t));
    -y.atan2(x)
}

/// Calculates the unsigned dihedral angle, in the range 0 to π, of the four ordered coordinates
///
/// # Examples
///
/// ```
/// use dihedral::dihedral_unsigned;
///
/// let P0 = [24.969, 13.428, 30.692]; // N
/// let P1 = [24.044, 12.661, 29.808]; // CA
/// let P2 = [22.785, 13.482, 29.543]; // C
/// let P3 = [21.951, 13.670, 30.431]; // O
/// let P4 = [23.672, 11.328, 30.466]; // CB
/// let P5 = [22.881, 10.326, 29.620]; // CG
/// let P6 = [23.691, 9.935, 28.389]; // CD1
/// let P7 = [22.557, 9.096, 30.459]; // CD2
///
/// assert!((dihedral_unsigned([P0, P1, P2, P3]).to_degrees() - (71.21515)).abs() < 1E-2);
/// assert!((dihedral_unsigned([P0, P1, P4, P5]).to_degrees() - (171.94319)).abs() < 1E-2);
/// assert!((dihedral_unsigned([P1, P4, P5, P6]).to_degrees() - (60.82226)).abs() < 1E-2);
/// assert!((dihedral_unsigned([P1, P4, P5, P7]).to_degrees() - (177.63641)).abs() < 1E-2);
/// ```
#[inline]
pub fn dihedral_unsigned([a, b, c, d]: [Point3D; 4]) -> Float {
    let (a, b, c) = (v(a, b), v(b, c), v(c, d));
    let (r, s) = (cross(a, b), cross(b, c));
    (dot(r, s) / (norm(r) * norm(s))).acos()
}

/// Norm (length) of a 3D vector
#[inline(always)]
fn norm(a: Vector3D) -> Float {
    dot(a, a).sqrt()
}

/// Unit vector
#[inline(always)]
fn u(v: Vector3D) -> Vector3D {
    let norm = norm(v);
    [v[0] / norm, v[1] / norm, v[2] / norm]
}

/// The vector connecting the two points
#[inline(always)]
fn v(p1: Point3D, p2: Point3D) -> Vector3D {
    [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]]
}

/// Dot product between two 3D vectors
#[inline(always)]
fn dot(m: Vector3D, n: Vector3D) -> Float {
    m[0] * n[0] + m[1] * n[1] + m[2] * n[2]
}

/// Cross product between two 3D vectors
#[inline(always)]
fn cross(m: Vector3D, n: Vector3D) -> Vector3D {
    [
        m[1] * n[2] - m[2] * n[1],
        m[2] * n[0] - m[0] * n[2],
        m[0] * n[1] - m[1] * n[0],
    ]
}

#[cfg(test)]
mod tests {

    use super::*;

    const P0: Point3D = [24.969, 13.428, 30.692]; // N
    const P1: Point3D = [24.044, 12.661, 29.808]; // CA
    const P2: Point3D = [22.785, 13.482, 29.543]; // C
    const P3: Point3D = [21.951, 13.670, 30.431]; // O
    const P4: Point3D = [23.672, 11.328, 30.466]; // CB
    const P5: Point3D = [22.881, 10.326, 29.620]; // CG
    const P6: Point3D = [23.691, 9.935, 28.389]; // CD1
    const P7: Point3D = [22.557, 9.096, 30.459]; // CD2

    #[test]
    fn test_dihedral() {
        assert!((dihedral([P0, P1, P2, P3]).to_degrees() - (-71.21515)).abs() < 1E-2);
        assert!((dihedral([P0, P1, P4, P5]).to_degrees() - (-171.94319)).abs() < 1E-2);
        assert!((dihedral([P1, P4, P5, P6]).to_degrees() - (60.82226)).abs() < 1E-2);
        assert!((dihedral([P1, P4, P5, P7]).to_degrees() - (-177.63641)).abs() < 1E-2);
    }

    #[test]
    fn test_dihedral_unsigned() {
        println!("{}", dihedral_unsigned([P0, P1, P4, P5]).to_degrees());
        assert!((dihedral_unsigned([P0, P1, P2, P3]).to_degrees() - (71.21515)).abs() < 1E-2);
        assert!((dihedral_unsigned([P0, P1, P4, P5]).to_degrees() - (171.94319)).abs() < 1E-2);
        assert!((dihedral_unsigned([P1, P4, P5, P6]).to_degrees() - (60.82226)).abs() < 1E-2);
        assert!((dihedral_unsigned([P1, P4, P5, P7]).to_degrees() - (177.63641)).abs() < 1E-2);
    }
}
