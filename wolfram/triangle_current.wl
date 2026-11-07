(* Independent planar surface-current Biot-Savart references.
   Run: wolframscript -file testing/wolfram/triangle_current.wl *)
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "common.wl"}]];
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "surface_current_common.wl"}]];

triangleVertices = {
  {-1/2, -1/2, -1/4}, {1/2, -1/2, -1/4}, {0, 1/2, -1/4}
};

triangleCases = {
  {"normal_above", {1, 0, 0}, {0, 0, 3/4}},
  {"reverse_transverse", {0, -3/2, 0}, {1/8, 0, 1/2}},
  {"mixed_offaxis", {2/3, -1/4, 3/5}, {2/5, 1/5, 3/4}},
  {"below", {0, 2, 0}, {0, -1/8, -1/2}},
  {"near_face_above", {1, 1/2, 0}, {0, -1/8, -1/4 + 1/64}},
  {"near_face_below", {-1/3, 2/5, 1/4}, {1/8, -1/16, -1/4 - 1/64}},
  {"near_edge", {1/2, -1/3, 0}, {0, -1/2 + 1/128, -1/4 + 1/64}},
  {"near_vertex", {3/2, 0, 0}, {-1/2 + 1/128, -1/2 + 1/128, -1/4 + 1/32}},
  {"pure_normal", {0, 0, 7/3}, {0, 0, 1/2}},
  {"far_field", {7/5, 1/3, 0}, {3, 2, 4}}
};

triangleRows = Table[
  Module[{label = item[[1]], k = item[[2]], p = item[[3]], field},
    field = hpConfirmed["triangle_current/" <> label,
      hpSurfaceCurrentField[k, p, triangleVertices, #] &];
    If[label == "pure_normal" && field != {0, 0, 0},
      Print["FAIL triangle_current normal density must project to zero"]; Exit[1]
    ];
    Join[{label}, hpASCII /@ k, hpASCII /@ p, hpASCII /@ field]
  ], {item, triangleCases}
];
Print["TriangleCurrent per-face tangent projection checked."];
hpWriteCSV["triangle_current",
  {"case", "kx", "ky", "kz", "x", "y", "z", "bx", "by", "bz"},
  triangleRows];
