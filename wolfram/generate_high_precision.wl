(* Independent magnetostatic reference values. Run with:
   wolframscript -file testing/wolfram/generate_high_precision.wl *)

ClearAll[ascii, confirmed, loopField, faceField, cuboidField, writeCases];

(* The exact decimal permeability used by Magba for current sources, H/m. *)
mu0 = 125663706127/10^17;
outputDirectory = FileNameJoin[{
   DirectoryName[ExpandFileName[$InputFileName]], "..", "data", "high-precision"
}];

(* RealDigits yields exactly 26 significant digits with no Mathematica-specific
   backtick precision syntax, including for arbitrary-precision zeros. *)
ascii[value_] := Module[{digits, exponent, sign},
  (* Residual symmetry components under 1e-40 T are numerical zero. *)
  If[PossibleZeroQ[value] || TrueQ[Abs[value] < 10^-40], Return["0"]];
  {digits, exponent} = RealDigits[Abs[value], 10, 26];
  If[!VectorQ[digits, IntegerQ],
    Print["FAIL cannot format a value with 26 significant digits: ", InputForm[value]];
    Exit[1]
  ];
  If[Total[digits] == 0, Return["0"]];
  sign = If[TrueQ[value < 0], "-", ""];
  sign <> ToString[First[digits]] <> "." <>
    StringJoin[ToString /@ Rest[digits]] <> "e" <> ToString[exponent - 1]
];

(* Independently integrate d\[ScriptL] x (r-r') / |r-r'|^3 along the wire. *)
loopField[diameter_, current_, p_, wp_] := Module[{radius = diameter/2},
  With[{a = radius, x = p[[1]], y = p[[2]], z = p[[3]]},
    mu0 current/(4 Pi) NIntegrate[
      Evaluate[a {z Cos[t], z Sin[t], a - x Cos[t] - y Sin[t]}/
        (a^2 + x^2 + y^2 + z^2 - 2 a (x Cos[t] + y Sin[t]))^(3/2)],
      {t, 0, 2 Pi}, WorkingPrecision -> wp,
      AccuracyGoal -> wp - 12, PrecisionGoal -> wp - 12,
      Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
      MaxRecursion -> 20
    ]
  ]
];

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

(* J is magnetic polarization in teslas, so mu0 cancels in the surface term.
   B = J + mu0 H inside, B = mu0 H outside; the boundary is excluded. *)
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

(* Two independent precision settings, with a componentwise absolute error
   normalized by the overall field magnitude to accommodate symmetry zeros. *)
confirmed[label_, evaluator_] := Module[{low, high, error, scale},
  low = Check[evaluator[40], $Failed];
  high = Check[evaluator[60], $Failed];
  If[!VectorQ[low, NumericQ] || !VectorQ[high, NumericQ] ||
     Length[low] != 3 || Length[high] != 3,
    Print["FAIL ", label, ": nonnumeric integral"]; Exit[1]
  ];
  error = Max[Abs[high - low]];
  scale = Max[10^-18, Max[Abs[high]]];
  Print[label, ": scaled 40/60 disagreement ",
    If[PossibleZeroQ[error/scale], "below working-precision resolution",
      ascii[error/scale]], " (limit 1e-22)"];
  If[!TrueQ[error < 10^-22 scale],
    Print["FAIL ", label, ": 40/60 precision disagreement"]; Exit[1]
  ];
  high
];

loopCases = {
  {"axis_center", 2, 3, {0, 0, 0}},
  {"axis_above", 2, 3, {0, 0, 1/2}},
  {"axis_below_reverse", 3, -2, {0, 0, -3/4}},
  {"radial_inside_plane", 2, 2, {1/2, 0, 0}},
  {"radial_outside_plane", 2, -3, {3/2, 0, 0}},
  {"offaxis_above", 2, 1, {3/5, 2/5, 7/10}},
  {"offaxis_below", 3, 5, {-1/4, 1/2, -3/8}},
  {"near_wire_inner", 2, 1, {127/128, 1/64, 1/128}},
  {"near_wire_outer", 2, 1, {129/128, -1/64, -1/128}},
  {"unequal_xy", 5/4, 7/3, {2/5, -3/10, 9/10}},
  {"oblique_reverse", 4/3, -3/2, {-7/10, 1/2, 1/5}},
  {"millimeter_scale", 1/100, 1/2, {1/1000, -1/500, 1/500}}
};

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

Print["Generating high-precision independent reference values with Wolfram Engine."];
Print["Exact mu0 = ", ascii[mu0], " H/m; WorkingPrecision = 40 and 60."];

loopRows = Table[
  Module[{label = item[[1]], d = item[[2]], i = item[[3]],
          p = item[[4]], b, exact},
    b = confirmed["loop/" <> label, loopField[d, i, p, #] &];
    If[p[[1]] == 0 && p[[2]] == 0,
      exact = mu0 i (d/2)^2/(2 ((d/2)^2 + p[[3]]^2)^(3/2));
      If[!TrueQ[Abs[b[[3]] - exact] < 10^-23 Max[10^-18, Abs[exact]]] ||
         !TrueQ[Max[Abs[b[[1 ;; 2]]]] < 10^-23 Max[10^-18, Abs[exact]]],
        Print["FAIL ", label, ": analytic on-axis Biot-Savart check"]; Exit[1]
      ];
      Print["  analytic on-axis Bz check passed"]
    ];
    Join[{label}, ascii /@ {d, i}, ascii /@ p, ascii /@ b]
  ], {item, loopCases}
];

cuboidRows = Table[
  Module[{label = item[[1]], j = item[[2]],
          d = item[[3]], p = item[[4]], b},
    b = confirmed["cuboid/" <> label, cuboidField[j, d, p, #] &];
    If[label == "cube_center_z",
      If[!TrueQ[Abs[b[[3]] - 2/3] < 10^-22],
        Print["FAIL cube-center demagnetizing factor: Bz = ", ascii[b[[3]]],
          " T != 2/3 T"]; Exit[1]
      ];
      Print["  cube-center Bz = 2/3 T check passed"]
    ];
    Join[{label}, ascii /@ j, ascii /@ d, ascii /@ p, ascii /@ b]
  ], {item, cuboidCases}
];

If[!DirectoryQ[outputDirectory],
  CreateDirectory[outputDirectory, CreateIntermediateDirectories -> True]
];
Export[FileNameJoin[{outputDirectory, "circular_current.csv"}],
  Prepend[loopRows, {"case", "diameter", "current", "x", "y", "z", "bx", "by", "bz"}],
  "CSV"];
Export[FileNameJoin[{outputDirectory, "cuboid.csv"}],
  Prepend[cuboidRows, {"case", "jx", "jy", "jz", "dx", "dy", "dz", "x", "y", "z", "bx", "by", "bz"}],
  "CSV"];
Print["Wrote ", Length[loopRows], " loop and ", Length[cuboidRows],
  " cuboid cases to ", outputDirectory];
