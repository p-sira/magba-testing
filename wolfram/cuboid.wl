(* Independent Wolfram reference generator for cuboid.
   Run: wolframscript -file wolfram/cuboid.wl *)

Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "common.wl"}]];

ClearAll[faceField, cuboidField];

(* The face's magnetic surface-charge field is a two-dimensional integral of
   (point-source)/|point-source|^3. Integrate source coordinate u analytically
   and use one-dimensional high-precision quadrature on v. This avoids the
   stagnating error estimates of generic adaptive 2D cubature. *)
faceField[p_, dimensions_, axis_, side_, wp_] := Module[
  {others = DeleteCases[Range[3], axis], primitive},
  With[{i = others[[1]], j = others[[2]], k = axis,
        c = side dimensions[[axis]]/2, d = dimensions, point = p},
    primitive[uValue_] := With[
      {t = uValue - point[[i]], q = (point[[j]] - v)^2 + (point[[k]] - c)^2},
      ReplacePart[{0, 0, 0}, {
        i -> 1/Sqrt[t^2 + q],
        j -> (point[[j]] - v) t/(q Sqrt[t^2 + q]),
        k -> (point[[k]] - c) t/(q Sqrt[t^2 + q])
      }]
    ];
    NIntegrate[
      Evaluate[primitive[d[[i]]/2] - primitive[-d[[i]]/2]],
      {v, -d[[j]]/2, d[[j]]/2},
      WorkingPrecision -> wp, AccuracyGoal -> wp - 12,
      PrecisionGoal -> wp - 12,
      Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
      MaxRecursion -> 20
    ]
  ]
];

(* J is magnetic polarization in teslas, so hpMu0 cancels in the surface term.
   B = J + hpMu0 H inside, B = hpMu0 H outside; the boundary is excluded. *)
cuboidField[j_, dimensions_, p_, wp_] := Module[
  {field = {0, 0, 0}, inside},
  Do[
    If[j[[axis]] != 0,
      Do[field += side j[[axis]] faceField[p, dimensions, axis, side, wp],
        {side, {-1, 1}}]
    ],
    {axis, 1, 3}
  ];
  inside = And @@ MapThread[Abs[#1] < #2/2 &, {p, dimensions}];
  (field/(4 Pi)) + If[inside, j, {0, 0, 0}]
];


cuboidCases = {
  {"cube_center_z", {0, 0, 1}, {2, 2, 2}, {0, 0, 0}},
  {"cube_center_vector", {1/3, -1/5, 2/3}, {2, 2, 2}, {0, 0, 0}},
  {"cube_axis_outside", {0, 0, 1}, {2, 2, 2}, {0, 0, 3}},
  {"cube_far_oblique", {0, 0, 1}, {2, 2, 2}, {3, 2, 5}},
  {"rect_center_z", {0, 0, 1}, {1, 2, 3}, {0, 0, 0}},
  {"rect_inside_oblique", {2/5, -1/3, 4/5}, {1, 2, 3}, {1/5, -1/4, 1/3}},
  {"rect_outside_side", {3/4, 0, 0}, {1, 2, 3}, {4/5, 1/5, 1/4}},
  {"rect_outside_corner", {1/3, -2/5, 3/4}, {1, 2, 3}, {4/5, 6/5, 19/10}},
  {"near_face_outer", {0, 0, 1}, {2, 2, 2}, {1/4, 1/8, 1025/1024}},
  {"near_face_inner", {0, 0, 1}, {2, 2, 2}, {1/4, 1/8, 1023/1024}},
  {"near_edge", {0, 0, 1}, {2, 2, 2}, {63/64, 63/64, 5/4}},
  {"negative_polarization", {1/3, 1/5, -1/2}, {3/2, 5/4, 2}, {-7/10, 1/8, -11/10}}
};


cuboidRows = Table[
  Module[{label = item[[1]], j = item[[2]],
          d = item[[3]], p = item[[4]], b},
    b = hpConfirmed["cuboid/" <> label, cuboidField[j, d, p, #] &];
    If[label == "cube_center_z",
      If[!TrueQ[Abs[b[[3]] - 2/3] < 10^-22],
        Print["FAIL cube-center demagnetizing factor: Bz = ", hpASCII[b[[3]]],
          " T != 2/3 T"]; Exit[1]
      ];
      Print["  cube-center Bz = 2/3 T check passed"]
    ];
    Join[{label}, hpASCII /@ j, hpASCII /@ d, hpASCII /@ p, hpASCII /@ b]
  ], {item, cuboidCases}
];


hpWriteCSV["cuboid", {"case", "jx", "jy", "jz", "dx", "dy", "dz", "x", "y", "z", "bx", "by", "bz"}, cuboidRows];
