(* Independent high-precision open TriangleMagnet surface-charge reference.
   Run: wolframscript -file testing/wolfram/triangle_magnet.wl *)
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "common.wl"}]];
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "polyhedron_common.wl"}]];

altTriangle = polyTriangleVerticesAlt;
altTriangleCenter = Total[altTriangle]/3;
altTriangleNormal = Cross[
  altTriangle[[2]] - altTriangle[[1]],
  altTriangle[[3]] - altTriangle[[1]]
];
altTriangleEdge = (altTriangle[[1]] + altTriangle[[2]])/2;

triangleCases = {
  {"above_centroid", 0, {0, 0, 1}, {-1/4, -1/6, 1/4}},
  {"below_centroid", 0, {0, 0, 1}, {1/4, -1/6, -1/4}},
  {"near_face_above", 0, {1/3, -1/5, 2/3}, {-1/128, -1/6, 1/128}},
  {"near_face_below", 0, {-1/2, 1/4, 1/3}, {1/128, -1/6, -1/128}},
  {"above_offaxis", 0, {1, 0, 0}, {1/8, -1/10, 1/2}},
  {"outside_edge", 0, {1/2, 1/3, -1/4}, {1/2, -3/4, 5/8}},
  {"outside_vertex", 0, {-3/5, 2/3, 3/7}, {-5/8, -5/8, -3/8}},
  {"far_oblique", 0, {1/3, -2/5, 4/5}, {-1, 3/4, 3/5}},
  {"shape1_near_face_positive", 1, {1/3, -2/5, 4/5},
    altTriangleCenter + altTriangleNormal/256},
  {"shape1_near_face_negative", 1, {-1/2, 1/3, -3/5},
    altTriangleCenter - altTriangleNormal/256},
  {"shape1_near_edge_inner", 1, {3/4, -1/3, 1/2},
    altTriangleEdge + (altTriangle[[3]] - altTriangleEdge)/64 +
      altTriangleNormal/128},
  {"shape1_near_edge_outer", 1, {-2, 1, 1/3},
    altTriangleEdge - (altTriangle[[3]] - altTriangleEdge)/64 -
      altTriangleNormal/64},
  {"shape1_far_weak", 1, {1/1024, -1/512, 1/256}, {3, -2, 5/2}}
};

If[!TrueQ[Cross[polyTriangleVertices[[2]] - polyTriangleVertices[[1]],
    polyTriangleVertices[[3]] - polyTriangleVertices[[1]]] == {-1, 0, 1}],
  Print["FAIL triangle orientation"]; Exit[1]
];
polyGenerate["triangle_magnet", triangleCases, polyTriangleGeometries,
  polyTriangleB[polyTriangleGeometries[[#1 + 1, 1]], #2, #3, #4] &];
