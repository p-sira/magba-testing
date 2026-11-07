# High-precision magnetic-field references

Run from the `magba-testing` repository root:

```sh
wolframscript -file wolfram/generate_high_precision.wl
```

The script deterministically writes `data/high-precision/circular_current.csv`
and `data/high-precision/cuboid.csv`. No Wolfram dependency is needed to read
the checked-in CSVs or run Magba's Rust tests.

Both geometries are centered at the origin with identity orientation, and all
coordinates and dimensions are in meters. The current loop lies in the XY
plane, has `diameter` in meters, and carries `current` in amperes (positive
current travels counterclockwise as viewed from +Z). Its magnetic flux density
uses the numerical Biot–Savart line integral and Magba's **exact** decimal
`mu0 = 1.25663706127e-6` H/m, rather than `4 Pi 10^-7`. Cuboids have side
lengths `dx,dy,dz` in meters and uniform magnetic polarization `jx,jy,jz`
in teslas. Their field independently integrates magnetic surface charge over
the six rectangular faces, adding the internal polarization `J` for points
strictly inside the cuboid: `B = J + mu0 H` inside and `B = mu0 H` outside.
For this formulation `mu0` cancels from the cuboid surface contribution. Each
face starts from the two-dimensional Coulomb surface integral; one source
coordinate is integrated analytically, and the remaining coordinate is
integrated numerically for reliable arbitrary precision.

All generator inputs are exact Wolfram integers or rationals, and outputs are
base-10 ASCII floating-point values with 26 significant digits that parse in
Rust. Numerical symmetry residuals below `1e-40` T are written as zero.
Each field is evaluated
by numerical quadrature at both 40 and 60 digits `WorkingPrecision`, using
`AccuracyGoal` and `PrecisionGoal` equal to `WorkingPrecision - 12`.
The maximum componentwise absolute difference must be less than `1e-22` times
the greater of the field's largest component magnitude and `1e-18` T.
Diagnostics display this normalized error for every case. The script also
checks the exact on-axis loop formula and the cube-center demagnetizing result
`Bz = 2/3 T` for `J = (0,0,1) T`.

Fixtures cover axes, transverse/off-axis points, reversed current, geometric
scales, interior and exterior cuboid points, multiple polarization components,
and near-wire/near-face cases. Most inputs are dyadic rationals, exactly
representable in IEEE `f64`; a few non-dyadic cases also test general inputs.
All samples avoid exact wires, faces, edges,
and corners, where singular integral evaluations or field conventions require
special treatment. The reference convergence threshold is a generator
consistency check; it is stricter than a typical `f64` comparison tolerance.
