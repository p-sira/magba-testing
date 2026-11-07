(* Independent finite-wire Biot-Savart references.
   Run: wolframscript -file testing/wolfram/path_current.wl *)
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "common.wl"}]];

ClearAll[pathSegmentField, pathField, pathSegmentExact];
pathVertices = {
  {
    {-1/2, -1/4, -1/4}, {1/2, -1/4, -1/4},
    {0, 1/2, -1/4}, {0, 0, 1/2}
  },
  {
    {-3/4, -1/2, 1/4}, {1/2, -1/4, 1/4},
    {3/4, 1/2, -1/4}, {-1/4, 3/4, 1/2}
  }
};

(* Direct Biot-Savart line integration, independent of finite-wire kernels. *)
pathSegmentField[a_, b_, p_, wp_] := With[{d = b - a, r = p - a},
  NIntegrate[
    Evaluate[Cross[d, r - t d]/((r - t d).(r - t d))^(3/2)],
    {t, 0, 1}, WorkingPrecision -> wp,
    AccuracyGoal -> wp - 12, PrecisionGoal -> wp - 12,
    Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
    MaxRecursion -> 24
  ]
];
pathField[current_, p_, vertices_, wp_] :=
  If[current == 0, {0, 0, 0},
    hpMu0 current/(4 Pi) Total[
      MapThread[pathSegmentField[#1, #2, p, wp] &,
        {Most[vertices], Rest[vertices]}]
    ]
  ];

(* Independent closed-form integral used only to check the quadrature. *)
pathSegmentExact[a_, b_, p_] := Module[
  {d = b - a, r1 = p - a, r2 = p - b, denom},
  denom = (d.d) (r1.r1) - (d.r1)^2;
  Cross[d, r1] ((d.r1)/Sqrt[r1.r1] - (d.r2)/Sqrt[r2.r2])/denom
];

pathCases = {
  {"far_oblique", 0, 3/2, {3/2, 5/4, 7/4}},
  {"near_first_wire", 0, 2, {0, -1/4 + 1/32, -1/4 + 1/64}},
  {"below_first_wire", 0, -3/4, {1/8, -3/8, -1/2}},
  {"near_second_wire", 0, 1, {1/4, 1/8, -1/4 + 1/32}},
  {"near_last_vertex", 0, 1, {1/128, 1/256, 1/2 + 1/64}},
  {"side_offaxis", 0, 5/4, {-3/8, 1/4, 1/8}},
  {"above_reverse", 0, -7/3, {0, 0, 3/4}},
  {"below", 0, 7/5, {0, 0, -3/4}},
  {"zero_current", 0, 0, {1/4, 1/4, 1/4}},
  {"shape1_near_first_wire", 1, 9/4, {-1/8, -3/8, 1/4 + 1/64}},
  {"shape1_near_last_wire", 1, -3/2, {1/4, 5/8, 1/8 + 1/64}},
  {"shape1_near_joint", 1, 3/4, {1/2 + 1/64, -1/4 + 1/128, 1/4 + 1/64}},
  {"shape1_far_scale", 1, 7/5, {7/2, -11/4, 9/4}}
};

pathRows = Table[
  Module[{label = item[[1]], shape = item[[2]], current = item[[3]],
          p = item[[4]], vertices, field, exact},
    vertices = pathVertices[[shape + 1]];
    field = hpConfirmed["path_current/" <> label,
      pathField[current, p, vertices, #] &];
    exact = hpMu0 current/(4 Pi) Total[
      MapThread[pathSegmentExact[#1, #2, p] &,
        {Most[vertices], Rest[vertices]}]
    ];
    If[!TrueQ[Max[Abs[field - N[exact, 60]]] <
        10^-22 Max[10^-18, Max[Abs[field]]]],
      Print["FAIL path_current/", label, ": analytic finite-wire check"]; Exit[1]
    ];
    Join[{label, ToString[shape]}, hpASCII /@ {current},
      hpASCII /@ p, hpASCII /@ field]
  ], {item, pathCases}
];
Print["PathCurrent analytic finite-wire checks passed."];
hpWriteCSV["path_current",
  {"case", "shape", "current", "x", "y", "z", "bx", "by", "bz"}, pathRows];
