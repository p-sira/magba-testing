(* Point-dipole references from the vacuum magnetostatic Green function.
   Moments are in A m^2, coordinates in m, and B in T. *)
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "common.wl"}]];

ClearAll[dipoleField];
dipoleField[m_, p_, precision_] := Module[{r2 = p.p, mr = m.p},
  If[r2 == 0, Print["FAIL dipole observer at origin"]; Exit[1]];
  N[hpMu0/(4 Pi) (3 p mr/r2^(5/2) - m/r2^(3/2)), precision]
];

(* Exact rational inputs; all sources are at the origin with no rotation. *)
dipoleCases = {
  {"axial_positive", {0, 0, 2}, {0, 0, 3/2}},
  {"axial_negative", {0, 0, -3}, {0, 0, -5/4}},
  {"axis_x", {5/4, 0, 0}, {2, 0, 0}},
  {"transverse_z", {0, 0, 1}, {2, 0, 0}},
  {"axis_x_reversed", {-1/3, 0, 0}, {-7/4, 0, 0}},
  {"transverse_y", {0, 5/6, 0}, {3/5, 0, 0}},
  {"offaxis_xy", {0, 5/6, 0}, {3/5, -7/10, 0}},
  {"offaxis_xyz", {2/3, -3/4, 5/6}, {7/5, -4/5, 3/2}},
  {"mixed_signs", {-2/5, 1/7, -5/3}, {-3/4, 1/2, -2/3}},
  {"cross_terms", {1/2, 1/3, 1/7}, {1/2, 2/3, 4/5}},
  {"near_origin", {1, -1, 2}, {1/100, 1/200, -1/100}},
  {"far_oblique", {5/2, -1/4, 1/3}, {15, -10, 7/2}},
  {"millimeter_scale", {1/1000000, -3/2000000, 1/3000000},
    {1/1000, -1/500, 3/1000}},
  {"zero_moment", {0, 0, 0}, {2/3, -3/4, 5/2}},
  {"micrometer_skew", {1/10^12, -2/10^12, 3/10^12},
    {1/10^6, -1/(2 10^6), 1/(4 10^6)}},
  {"micrometer_transverse", {-3/10^11, 1/10^11, 2/10^11},
    {0, 1/10^5, 0}},
  {"megameter_oblique", {3 10^6, -5 10^6, 7 10^6},
    {2 10^6, -10^6, 3 10^6}},
  {"submillimeter_oblique", {3/10000, -7/10000, 1/20000},
    {1/2000, -3/4000, 1/1000}}
};

dipoleRows = Table[
  Module[{label = item[[1]], m = item[[2]], p = item[[3]], b, expected},
    b = hpConfirmed["dipole/" <> label, dipoleField[m, p, #] &];
    If[StringStartsQ[label, "axial_"],
      expected = N[{0, 0, hpMu0 m[[3]]/(2 Pi Abs[p[[3]]]^3)}, 60];
      If[!TrueQ[Max[Abs[b - expected]] < 10^-45],
        Print["FAIL ", label, ": axial dipole identity"]; Exit[1]
      ]
    ];
    If[label == "transverse_z",
      expected = N[{0, 0, -hpMu0 m[[3]]/(4 Pi Abs[p[[1]]]^3)}, 60];
      If[!TrueQ[Max[Abs[b - expected]] < 10^-45],
        Print["FAIL ", label, ": transverse dipole identity"]; Exit[1]
      ]
    ];
    Join[{label}, hpASCII /@ m, hpASCII /@ p, hpASCII /@ b]
  ], {item, dipoleCases}
];

hpWriteCSV["dipole", {"case", "mx", "my", "mz", "x", "y", "z",
  "bx", "by", "bz"}, dipoleRows];
