(* Independent high-precision finite cylinder magnetic surface-charge fields.
   Run: wolframscript -file wolfram/cylinder.wl *)

Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "common.wl"}]];
ClearAll[hpCylinderDisk, hpCylinderSide, hpCylinderField];

(* The disk surface charge is +/-Jz. Integrate the radial source coordinate
   analytically in the Coulomb surface integral, leaving an angular quadrature.
   At fixed theta, a and b are the radial and tangential projections of (x,y).
   The primitives below integrate rho (p-r')/|p-r'|^3 over 0 <= rho <= R. *)
hpCylinderDisk[point_, radius_, faceZ_, wp_] := Module[
  {a, b, c, d2, radialPrimitive, scalarPrimitive, ir, is, expression},
  With[{x = point[[1]], y = point[[2]], z = point[[3]],
        r = radius, face = faceZ},
    a = x Cos[theta] + y Sin[theta];
    b = -x Sin[theta] + y Cos[theta];
    c = z - face;
    d2 = b^2 + c^2;
    radialPrimitive[t_] := -ArcSinh[t/Sqrt[d2]] + (t + a)/Sqrt[t^2 + d2];
    scalarPrimitive[t_] := -1/Sqrt[t^2 + d2] +
      a t/(d2 Sqrt[t^2 + d2]);
    ir = radialPrimitive[r - a] - radialPrimitive[-a];
    is = scalarPrimitive[r - a] - scalarPrimitive[-a];
    expression = {ir Cos[theta] - is b Sin[theta],
      ir Sin[theta] + is b Cos[theta], is c};
    NIntegrate[Evaluate[expression], {theta, 0, 2 Pi},
      WorkingPrecision -> wp, AccuracyGoal -> wp - 12,
      PrecisionGoal -> wp - 12,
      Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
      MaxRecursion -> 25]
  ]
];

(* Side charge is Jx cos(theta)+Jy sin(theta). Integrate the axial source
   coordinate z' analytically from -height/2 to height/2, leaving theta only. *)
hpCylinderSide[j_, point_, radius_, height_, wp_] := Module[
  {dx, dy, q, lower, upper, transverse, axial, expression},
  With[{x = point[[1]], y = point[[2]], z = point[[3]],
        r = radius, h = height, jx = j[[1]], jy = j[[2]]},
    dx = x - r Cos[theta];
    dy = y - r Sin[theta];
    q = dx^2 + dy^2;
    lower = -h/2 - z;
    upper = h/2 - z;
    transverse = upper/(q Sqrt[q + upper^2]) -
      lower/(q Sqrt[q + lower^2]);
    axial = 1/Sqrt[q + upper^2] - 1/Sqrt[q + lower^2];
    expression = r (jx Cos[theta] + jy Sin[theta]) *
      {dx transverse, dy transverse, axial};
    NIntegrate[Evaluate[expression], {theta, 0, 2 Pi},
      WorkingPrecision -> wp, AccuracyGoal -> wp - 12,
      PrecisionGoal -> wp - 12,
      Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
      MaxRecursion -> 25]
  ]
];

(* B = J + surface-charge integral/(4 Pi) inside, and just the surface
   integral/(4 Pi) outside. J is polarization in teslas, so mu0 cancels. *)
hpCylinderField[j_, diameter_, height_, point_, wp_] := Module[
  {r = diameter/2, total = {0, 0, 0}, inside},
  If[j[[3]] != 0,
    total += j[[3]] (hpCylinderDisk[point, r, height/2, wp] -
      hpCylinderDisk[point, r, -height/2, wp])
  ];
  If[j[[1]] != 0 || j[[2]] != 0,
    total += hpCylinderSide[j, point, r, height, wp]
  ];
  inside = point[[1]]^2 + point[[2]]^2 < r^2 && Abs[point[[3]]] < height/2;
  total/(4 Pi) + If[inside, j, {0, 0, 0}]
];

(* Uniform polarization cases with rational, usually dyadic, inputs. All
   coordinates are strictly away from the material boundary. *)
cylinderCases = {
  {"axial_center", {0, 0, 1}, 2, 2, {0, 0, 0}},
  {"axial_inside_axis", {0, 0, 3/4}, 2, 2, {0, 0, 3/4}},
  {"axial_outside_axis", {0, 0, 1}, 2, 2, {0, 0, 2}},
  {"axial_offaxis_inside", {0, 0, 1}, 2, 2, {1/4, -1/8, 3/8}},
  {"axial_near_axis", {0, 0, -1}, 2, 2, {1/1024, -1/512, 3/4}},
  {"axial_near_side_outer", {0, 0, 1}, 2, 2, {129/128, 1/64, 1/4}},
  {"axial_near_side_inner", {0, 0, 1}, 2, 2, {127/128, 1/64, -1/4}},
  {"diametric_center_x", {1, 0, 0}, 2, 2, {0, 0, 0}},
  {"diametric_center_y", {0, -1/2, 0}, 2, 2, {0, 0, 0}},
  {"diametric_offaxis_outside", {3/4, 1/4, 0}, 2, 2, {2, 1/2, 1/4}},
  {"mixed_inside", {1/2, -1/4, 3/4}, 3/2, 5/4, {1/4, -1/8, 1/8}},
  {"mixed_far", {-1/3, 2/5, -1/2}, 3/2, 3, {2, -2, 4}},
  (* Thin disk: diameter/height = 16, one point near its top face
     from inside and one outside its cylindrical rim. *)
  {"thin_disk_axial_near_face_inside", {0, 0, 1}, 2, 1/8, {1/4, 0, 7/128}},
  {"thin_disk_mixed_near_rim_outside", {1/2, -1/4, 3/4}, 2, 1/8,
    {129/128, 0, 1/32}},
  (* Long rod: height/diameter = 16, with near-side interior and
     near-cap exterior samples. *)
  {"long_rod_diametric_near_side_inside", {1, 0, 0}, 1/4, 4,
    {63/512, 0, 3/4}},
  {"long_rod_axial_near_cap_outside", {0, 0, -3/4}, 1/4, 4,
    {1/32, 1/64, 2 + 1/512}},
  (* A smaller thin disk exercises the side charge near the rim. *)
  {"small_disk_diametric_near_rim_outside", {-1/2, 0, 0}, 1/4, 1/64,
    {65/512, 0, 0}}
};

Print["Cylinder reference: exact rational geometry and polarization (T); ",
  "40/60-digit disk and side surface-charge quadrature."];
rows = Table[
  Module[{label = item[[1]], j = item[[2]], diameter = item[[3]],
          height = item[[4]], point = item[[5]], b, onAxis, radialCenter},
    b = hpConfirmed["cylinder/" <> label,
      hpCylinderField[j, diameter, height, point, #] &];
    If[point[[1]] == 0 && point[[2]] == 0 && j[[1]] == 0 && j[[2]] == 0,
      onAxis = j[[3]]/2 (
        (point[[3]] + height/2)/Sqrt[(diameter/2)^2 + (point[[3]] + height/2)^2] -
        (point[[3]] - height/2)/Sqrt[(diameter/2)^2 + (point[[3]] - height/2)^2]);
      If[!TrueQ[Max[Abs[b - {0, 0, onAxis}]] < 10^-22 Max[10^-18, Abs[onAxis]]],
        Print["FAIL ", label, ": analytic axial cylinder field"]; Exit[1]
      ];
      Print["  analytic axial cylinder Bz check passed"]
    ];
    If[point == {0, 0, 0} && j[[3]] == 0,
      (* Rotational symmetry and Nx+Ny+Nz=1 at an interior point give
         Bx/Jx = By/Jy = 1 - height/(4 Sqrt[R^2+(height/2)^2]). *)
      radialCenter = (1 - height/(4 Sqrt[(diameter/2)^2 + (height/2)^2])) j;
      If[!TrueQ[Max[Abs[b - radialCenter]] < 10^-22 Max[10^-18, Max[Abs[j]]]],
        Print["FAIL ", label, ": analytic diametric cylinder center field"]; Exit[1]
      ];
      Print["  analytic diametric center field check passed"]
    ];
    Join[{label}, hpASCII /@ j, hpASCII /@ {diameter, height},
      hpASCII /@ point, hpASCII /@ b]
  ], {item, cylinderCases}
];

hpWriteCSV["cylinder",
  {"case", "jx", "jy", "jz", "diameter", "height", "x", "y", "z", "bx", "by", "bz"},
  rows];
