use dihedral::dihedral;

type Point = [f64; 3];

const P0: Point = [24.969, 13.428, 30.692];
const P1: Point = [24.044, 12.661, 29.808];
const P2: Point = [22.785, 13.482, 29.543];
const P3: Point = [21.951, 13.670, 30.431];
const P4: Point = [23.672, 11.328, 30.466];
const P5: Point = [22.881, 10.326, 29.620];
const P6: Point = [23.691, 9.935, 28.389];
const P7: Point = [22.557, 9.096, 30.459];

fn main() {
    println!("{}", dihedral(&[P0, P1, P2, P3]));
    let seq = vec![P0, P1, P2, P3, P4, P5, P6, P7];
    seq.windows(4)
        .map(|frame| println!("{}", dihedral(&[frame[0], frame[1], frame[2], frame[3]])))
        .for_each(drop);
}
