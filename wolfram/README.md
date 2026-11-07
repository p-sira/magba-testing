# Independent Wolfram field references

Generate any fixture independently from the `magba-testing` root:

```sh
wolframscript -file wolfram/cylinder.wl
```

Each `wolfram/<source>.wl` creates only its matching
`data/high-precision/<source>.csv`. Supported source names:
`circular_current`, `path_current`, `triangle_current`, `sheet_current`,
`dipole`, `sphere`, `cuboid`, `cylinder`, `triangle_magnet`,
`tetrahedron_magnet` and `mesh_magnet`.

Source geometry and test cases are declared with exact rational inputs in
each script. All sources have origin position and identity orientation.
Lengths are meters; current is amperes; surface current density is A/m;
dipole moment is A·m²; magnet polarization and computed B are teslas.
Wolfram is **not** required to run the Rust tests consuming the committed CSVs.

## Independent reference methods

- **Currents:** Biot–Savart wire and surface integrals; triangle/sheet
  current density is projected onto each face.
- **Dipole and sphere:** magnetostatic Green function, interior sphere
  demagnetizing field and exterior multipole potential.
- **Solid magnets:** magnetic surface-charge Coulomb integrals with an
  interior polarization term; triangular faces have no interior term.
  Some surface coordinates are integrated analytically before independent
  numerical quadrature. Tetrahedron and mesh use identical tetrahedral
  geometry to test distinct Magba implementations.

`common.wl` supplies the exact decimal permeability used by Magba,
`μ₀ = 1.25663706127 × 10⁻⁶` H/m, plus numeric checks and CSV formatting;
`polyhedron_common.wl` shares triangular surface integration methods.
Uniform-magnet polarization fields have no dependence on `μ₀` after
conversion to magnetic flux density.

Integrals are evaluated at 40 and 60 digits of working precision; the
componentwise difference must be below
`1e-20 × max(1e-18 T, largest |B component|)`. Individual generators check
additional analytical identities where available. Output is deterministic
ASCII decimal CSV with 26 significant digits; numerical symmetry residuals
below `1e-40 T` become zero. Regular interior/exterior, off-axis, different
scale and near-boundary cases are included; **exact singular boundaries are
excluded**.
