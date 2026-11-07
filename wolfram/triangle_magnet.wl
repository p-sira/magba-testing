(* Independent high-precision open TriangleMagnet surface-charge reference.
   Run: wolframscript -file testing/wolfram/triangle_magnet.wl *)
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "common.wl"}]];
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "polyhedron_common.wl"}]];

triangleCases = {
  {"above_centroid", {0, 0, 1}, {-1/4, -1/6, 1/4}},
  {"below_centroid", {0, 0, 1}, {1/4, -1/6, -1/4}},
  {"near_face_above", {1/3, -1/5, 2/3}, {-1/128, -1/6, 1/128}},
  {"near_face_below", {-1/2, 1/4, 1/3}, {1/128, -1/6, -1/128}},
  {"above_offaxis", {1, 0, 0}, {1/8, -1/10, 1/2}},
  {"outside_edge", {1/2, 1/3, -1/4}, {1/2, -3/4, 5/8}},
  {"outside_vertex", {-3/5, 2/3, 3/7}, {-5/8, -5/8, -3/8}},
  {"far_oblique", {1/3, -2/5, 4/5}, {-1, 3/4, 3/5}}
};

If[!TrueQ[Cross[polyTriangleVertices[[2]] - polyTriangleVertices[[1]],
    polyTriangleVertices[[3]] - polyTriangleVertices[[1]]] == {-1, 0, 1}],
  Print["FAIL triangle orientation"]; Exit[1]
];
polyGenerate["triangle_magnet", triangleCases, polyTriangleB];
