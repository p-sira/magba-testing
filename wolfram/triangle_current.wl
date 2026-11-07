(* Independent planar surface-current Biot-Savart references.
   Run: wolframscript -file testing/wolfram/triangle_current.wl *)
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "common.wl"}]];
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "surface_current_common.wl"}]];

triangleVertices = {
  {
    {-1/2, -1/2, -1/4}, {1/2, -1/2, -1/4}, {0, 1/2, -1/4}
  },
  {
    {-3/4, -1/4, 1/3}, {1/2, -1/2, 1/8}, {1/4, 3/4, -1/3}
  }
};
triangleNormal1 = Cross[
  triangleVertices[[2, 2]] - triangleVertices[[2, 1]],
  triangleVertices[[2, 3]] - triangleVertices[[2, 1]]
];
triangleCentroid1 = Total[triangleVertices[[2]]]/3;
triangleEdgeMidpoint1 = (triangleVertices[[2, 1]] + triangleVertices[[2, 2]])/2;

triangleCases = {
  {"normal_above", 0, {1, 0, 0}, {0, 0, 3/4}},
  {"reverse_transverse", 0, {0, -3/2, 0}, {1/8, 0, 1/2}},
  {"mixed_offaxis", 0, {2/3, -1/4, 3/5}, {2/5, 1/5, 3/4}},
  {"below", 0, {0, 2, 0}, {0, -1/8, -1/2}},
  {"near_face_above", 0, {1, 1/2, 0}, {0, -1/8, -1/4 + 1/64}},
  {"near_face_below", 0, {-1/3, 2/5, 1/4}, {1/8, -1/16, -1/4 - 1/64}},
  {"near_edge", 0, {1/2, -1/3, 0}, {0, -1/2 + 1/128, -1/4 + 1/64}},
  {"near_vertex", 0, {3/2, 0, 0}, {-1/2 + 1/128, -1/2 + 1/128, -1/4 + 1/32}},
  {"pure_normal", 0, {0, 0, 7/3}, {0, 0, 1/2}},
  {"far_field", 0, {7/5, 1/3, 0}, {3, 2, 4}},
  {"shape1_near_face_plus", 1, {2/3, -1/3, 1/5},
    triangleCentroid1 + triangleNormal1/64},
  {"shape1_near_face_minus", 1, {-2/5, 1/2, 1/3},
    triangleCentroid1 - triangleNormal1/64},
  {"shape1_near_edge", 1, {3/4, 1/4, -1/2},
    triangleEdgeMidpoint1 + triangleNormal1/48},
  {"shape1_far_scale", 1, {7/3, -5/4, 1/6}, {7/2, -5/2, 11/3}}
};

triangleRows = Table[
  Module[{label = item[[1]], shape = item[[2]], k = item[[3]],
          p = item[[4]], field},
    field = hpConfirmed["triangle_current/" <> label,
      hpSurfaceCurrentField[k, p, triangleVertices[[shape + 1]], #] &];
    If[label == "pure_normal" && field != {0, 0, 0},
      Print["FAIL triangle_current normal density must project to zero"]; Exit[1]
    ];
    Join[{label, ToString[shape]}, hpASCII /@ k,
      hpASCII /@ p, hpASCII /@ field]
  ], {item, triangleCases}
];
Print["TriangleCurrent per-face tangent projection checked."];
hpWriteCSV["triangle_current",
  {"case", "shape", "kx", "ky", "kz", "x", "y", "z", "bx", "by", "bz"},
  triangleRows];
