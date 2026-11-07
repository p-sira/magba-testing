(* Uniformly polarized solid sphere: interior demagnetizing field and
   exterior scalar potential from the l=1 surface-charge multipole.
   Polarization J and B are in T, all lengths in m. *)
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "common.wl"}]];

ClearAll[sphereField];
sphereField[j_, diameter_, p_, precision_] := Module[
  {xx, yy, zz, coords, a = diameter/2, potential},
  If[diameter <= 0 || p.p == a^2,
    Print["FAIL sphere diameter or observer on surface"]; Exit[1]
  ];
  (* B = J + mu0 H = (2/3) J within a uniformly magnetized sphere. *)
  If[p.p < a^2, Return[N[2 j/3, precision]]];
  (* Outside, phi_m = a^3 (J.r)/(3 mu0 r^3). Taking -mu0 grad
     cancels mu0 and provides B independently of Magba's kernel. *)
  coords = {xx, yy, zz};
  potential = (a^3/3) (j.coords)/(coords.coords)^(3/2);
  N[-(D[potential, #] & /@ coords) /. Thread[coords -> p], precision]
];

(* Exact rational inputs, origin-centered spheres with identity orientation. *)
sphereCases = {
  {"center_z", {0, 0, 1}, 2, {0, 0, 0}},
  {"center_vector", {3/5, -2/7, 4/9}, 2, {0, 0, 0}},
  {"inside_oblique", {1/3, 2/5, -3/4}, 2, {1/4, -1/3, 2/5}},
  {"near_surface_inside", {2/3, -1/5, 1/2}, 2, {1023/1024, 0, 0}},
  {"axis_outside", {0, 0, 1}, 2, {0, 0, 2}},
  {"equator_outside", {0, 0, 1}, 2, {3/2, 0, 0}},
  {"outside_oblique", {1/3, -2/5, 3/4}, 2, {3/2, 2/3, -1/2}},
  {"near_surface_outside", {1/3, 1/2, -1/4}, 2, {1025/1024, 0, 0}},
  {"inside_mixed", {3/7, -2/5, 1/3}, 2, {1/2, 1/3, 1/4}},
  {"small_outside", {0, 2/3, 1/7}, 1/100, {1/100, 0, 1/200}},
  {"large_inside", {1/3, 4/5, -1/4}, 10, {-3, 2, 1}},
  {"far_oblique", {-1/5, 1/3, 2/7}, 3/2, {8, -5, 4}},
  {"tiny_near_surface_outside", {1/2, -3/5, 1/7}, 1/500,
    {1001/1000000, 0, 0}},
  {"zero_polarization", {0, 0, 0}, 2, {0, 0, 3/2}}
};

sphereRows = Table[
  Module[{label = item[[1]], j = item[[2]], d = item[[3]],
          p = item[[4]], b, expected},
    b = hpConfirmed["sphere/" <> label, sphereField[j, d, p, #] &];
    If[p.p < (d/2)^2,
      expected = N[2 j/3, 60];
      If[!TrueQ[Max[Abs[b - expected]] < 10^-45],
        Print["FAIL ", label, ": internal demagnetizing identity"]; Exit[1]
      ],
      (* External moment m = 4 Pi a^3 J/(3 mu0): verify the dipole limit. *)
      expected = N[((d/2)^3/3) (3 p (j.p)/(p.p)^(5/2) -
        j/(p.p)^(3/2)), 60];
      If[!TrueQ[Max[Abs[b - expected]] < 10^-45],
        Print["FAIL ", label, ": exterior dipole identity"]; Exit[1]
      ]
    ];
    Join[{label}, hpASCII /@ j, {hpASCII[d]}, hpASCII /@ p,
      hpASCII /@ b]
  ], {item, sphereCases}
];

hpWriteCSV["sphere", {"case", "jx", "jy", "jz", "diameter", "x",
  "y", "z", "bx", "by", "bz"}, sphereRows];
