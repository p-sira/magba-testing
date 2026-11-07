(* Coulomb surface-charge quadrature for planar triangles. All geometry and
   polarization inputs are exact rationals; no MagPyLib kernels are used.
   Integrate radially from vertex a analytically, leaving one numerical
   integration over rays to the opposite edge. *)

ClearAll[polyRayPrimitive, polyFaceKernel, polyTriangleB, polyClosedB,
  polyInsideQ, polyGenerate, polyTriangleVertices, polyTriangleVerticesAlt,
  polyTetraVertices, polyTetraVerticesAlt, polyTetraFaces,
  polyBipyramidVertices, polyBipyramidFaces, polyTriangleGeometries,
  polyTetraGeometries, polyMeshGeometries, polyOutwardQ,
  polyOffFacePlanesQ];

polyTriangleVertices = {
  {-1/2, -1/2, -1/2}, {1/2, -1/2, 1/2}, {0, 1/2, 0}
};
polyTriangleVerticesAlt = {
  {-3/4, -1/4, 1/3}, {1/2, -1/2, 1/8}, {1/4, 3/4, -1/3}
};
polyTetraVertices = {
  {-1/2, -1/2, -1/2}, {1/2, -1/2, -1/2},
  {0, 1/2, -1/2}, {0, 0, 1/2}
};
polyTetraVerticesAlt = {
  {-3/4, -1/2, -1/3}, {1/2, -1/4, -1/2},
  {1/4, 3/4, -1/4}, {-1/8, 1/8, 4/5}
};
(* Zero-based, outward-facing vertex indices as passed to Magba. *)
polyTetraFaces = {{0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {0, 3, 2}};
(* A distinct closed, five-vertex triangular bipyramid mesh. *)
polyBipyramidVertices = {
  {-3/4, -1/2, 0}, {3/4, -1/2, 0}, {0, 3/4, 0},
  {0, 0, 5/4}, {0, 0, -3/4}
};
polyBipyramidFaces = {
  {0, 1, 3}, {1, 2, 3}, {0, 3, 2},
  {0, 4, 1}, {1, 4, 2}, {0, 2, 4}
};
(* Shape IDs correspond exactly to the integer CSV column. *)
polyTriangleGeometries = {
  {polyTriangleVertices, {{0, 1, 2}}},
  {polyTriangleVerticesAlt, {{0, 1, 2}}}
};
polyTetraGeometries = {
  {polyTetraVertices, polyTetraFaces},
  {polyTetraVerticesAlt, polyTetraFaces}
};
polyMeshGeometries = Append[polyTetraGeometries,
  {polyBipyramidVertices, polyBipyramidFaces}];

(* Primitive of t (r-t d)/|r-t d|^3, at fixed ray direction d.
   h2 is positive provided the observer is off the face plane. *)
polyRayPrimitive[r_, d_, t_] := Module[
  {dd = d.d, rd = r.d, rr = r.r, along, h2, perpendicular, w, distance},
  along = rd/dd;
  h2 = rr - rd^2/dd;
  perpendicular = r - along d;
  w = t - along;
  distance = Sqrt[dd w^2 + h2];
  perpendicular (along w/h2 - 1/dd)/distance -
    d (ArcSinh[Sqrt[dd] w/Sqrt[h2]]/(dd Sqrt[dd]) -
       (w + along)/(dd distance))
];

(* s(u,t)=a+t ((1-u)(b-a)+u(c-a)); ds=|cross(b-a,c-a)| t du dt.
   The cross product and J dot n are applied outside this vector kernel. *)
polyFaceKernel[vertices_, p_, wp_] := Module[
  {a = vertices[[1]], e1, e2, r, u, ray},
  e1 = vertices[[2]] - a;
  e2 = vertices[[3]] - a;
  r = p - a;
  ray = (1 - u) e1 + u e2;
  NIntegrate[
    Evaluate[polyRayPrimitive[r, ray, 1] - polyRayPrimitive[r, ray, 0]],
    {u, 0, 1},
    WorkingPrecision -> wp, AccuracyGoal -> wp - 12,
    PrecisionGoal -> wp - 12,
    Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
    MaxRecursion -> 24
  ]
];

polyTriangleB[vertices_, j_, p_, wp_] := Module[
  {v = vertices, normalArea},
  normalArea = Cross[v[[2]] - v[[1]], v[[3]] - v[[1]]];
  (j.normalArea) polyFaceKernel[v, p, wp]/(4 Pi)
];

(* Every face must be outward-oriented. This strict test excludes points on
   the surface, for which the internal-field convention is ambiguous. *)
polyInsideQ[vertices_, faces_, p_] := And @@ Table[
  With[{v = vertices[[face + 1]]},
    Cross[v[[2]] - v[[1]], v[[3]] - v[[1]]].(p - v[[1]]) < 0
  ], {face, faces}
];

polyOutwardQ[vertices_, faces_] := Module[
  {center = Total[vertices]/Length[vertices]},
  And @@ Table[
    With[{v = vertices[[face + 1]]},
      Cross[v[[2]] - v[[1]], v[[3]] - v[[1]]].(center - v[[1]]) < 0
    ], {face, faces}
  ]
];

polyOffFacePlanesQ[vertices_, faces_, p_] := And @@ Table[
  With[{v = vertices[[face + 1]]},
    Cross[v[[2]] - v[[1]], v[[3]] - v[[1]]].(p - v[[1]]) != 0
  ], {face, faces}
];

polyClosedB[vertices_, faces_, j_, p_, wp_] := Module[
  {surface = {0, 0, 0}, v, normalArea},
  Do[
    v = vertices[[face + 1]];
    normalArea = Cross[v[[2]] - v[[1]], v[[3]] - v[[1]]];
    surface += (j.normalArea) polyFaceKernel[v, p, wp],
    {face, faces}
  ];
  surface/(4 Pi) + If[polyInsideQ[vertices, faces, p], j, {0, 0, 0}]
];

polyGenerate[stem_String, cases_List, geometries_List, field_] := Module[
  {rows, shape, vertices, faces, j, p},
  If[!DuplicateFreeQ[cases[[All, 1]]],
    Print["FAIL duplicate ", stem, " case labels"]; Exit[1]
  ];
  Do[
    shape = case[[2]]; j = case[[3]]; p = case[[4]];
    If[!IntegerQ[shape] || shape < 0 || shape >= Length[geometries],
      Print["FAIL invalid shape for ", stem, "/", case[[1]]]; Exit[1]
    ];
    {vertices, faces} = geometries[[shape + 1]];
    If[!TrueQ[polyOffFacePlanesQ[vertices, faces, p]],
      Print["FAIL on face plane for ", stem, "/", case[[1]]]; Exit[1]
    ];
    If[And @@ Table[
        j.Cross[vertices[[face[[2]] + 1]] - vertices[[face[[1]] + 1]],
          vertices[[face[[3]] + 1]] - vertices[[face[[1]] + 1]]] == 0,
        {face, faces}
      ],
      Print["FAIL zero magnetic surface charge for ", stem, "/", case[[1]]];
      Exit[1]
    ],
    {case, cases}
  ];
  Print["Generating ", stem, " from Coulomb surface charge (40/60 digits)"];
  rows = Table[
    With[{label = case[[1]], shape = case[[2]],
          j = case[[3]], p = case[[4]]},
      Join[{label, shape}, hpASCII /@ j, hpASCII /@ p,
        hpASCII /@ hpConfirmed[
          stem <> "/" <> label, field[shape, j, p, #] &]]
    ],
    {case, cases}
  ];
  hpWriteCSV[stem, {"case", "shape", "jx", "jy", "jz", "x", "y", "z",
     "bx", "by", "bz"}, rows]
];
