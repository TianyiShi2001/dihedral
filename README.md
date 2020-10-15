# dihedral

Dihedral angle calculations

This crate provides functions for working with dihedral angles. Currently, there are two functions:

- `dihedral` calculates the dihedral angle in the range -π to π in accordance with biochemistry textbooks (see also: https://en.wikipedia.org/wiki/Dihedral_angle#In_stereochemistry)
- `dihedral_unsigned` ignores the direction of rotation and outputs the angle within the range 0 to π. This function is faster than the above signed version.