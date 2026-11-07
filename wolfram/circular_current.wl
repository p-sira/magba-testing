(* Independent Wolfram reference generator for circular_current.
   Run: wolframscript -file wolfram/circular_current.wl *)

Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "common.wl"}]];

ClearAll[loopField];

(* Independently integrate d\[ScriptL] x (r-r') / |r-r'|^3 along the wire. *)
loopField[diameter_, current_, p_, wp_] := Module[{radius = diameter/2},
  With[{a = radius, x = p[[1]], y = p[[2]], z = p[[3]]},
    hpMu0 current/(4 Pi) NIntegrate[
      Evaluate[a {z Cos[t], z Sin[t], a - x Cos[t] - y Sin[t]}/
        (a^2 + x^2 + y^2 + z^2 - 2 a (x Cos[t] + y Sin[t]))^(3/2)],
      {t, 0, 2 Pi}, WorkingPrecision -> wp,
      AccuracyGoal -> wp - 12, PrecisionGoal -> wp - 12,
      Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
      MaxRecursion -> 20
    ]
  ]
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


loopRows = Table[
  Module[{label = item[[1]], d = item[[2]], i = item[[3]],
          p = item[[4]], b, exact},
    b = hpConfirmed["loop/" <> label, loopField[d, i, p, #] &];
    If[p[[1]] == 0 && p[[2]] == 0,
      exact = hpMu0 i (d/2)^2/(2 ((d/2)^2 + p[[3]]^2)^(3/2));
      If[!TrueQ[Abs[b[[3]] - exact] < 10^-23 Max[10^-18, Abs[exact]]] ||
         !TrueQ[Max[Abs[b[[1 ;; 2]]]] < 10^-23 Max[10^-18, Abs[exact]]],
        Print["FAIL ", label, ": analytic on-axis Biot-Savart check"]; Exit[1]
      ];
      Print["  analytic on-axis Bz check passed"]
    ];
    Join[{label}, hpASCII /@ {d, i}, hpASCII /@ p, hpASCII /@ b]
  ], {item, loopCases}
];


hpWriteCSV["circular_current", {"case", "diameter", "current", "x", "y", "z", "bx", "by", "bz"}, loopRows];
