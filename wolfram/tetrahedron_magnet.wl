(* Independent high-precision closed TetrahedronMagnet surface-charge reference.
   Run: wolframscript -file testing/wolfram/tetrahedron_magnet.wl *)
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "common.wl"}]];
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "polyhedron_common.wl"}]];

altTetra = polyTetraVerticesAlt;
altTetraCenter = Total[altTetra]/4;
altTetraBaseCenter = Total[altTetra[[polyTetraFaces[[1]] + 1]]]/3;
altTetraEdge = (altTetra[[1]] + altTetra[[2]])/2;

tetraCases = {
  {"interior_centroid", 0, {0, 0, 1}, {0, -1/8, -1/4}},
  {"interior_vector", 0, {1/3, -2/5, 3/4}, {1/16, -1/8, -1/4}},
  {"near_base_inside", 0, {0, 0, 1}, {0, -1/6, -63/128}},
  {"near_base_outside", 0, {0, 0, 1}, {0, -1/6, -65/128}},
  {"near_side_inside", 0, {-2/5, 1/3, 3/4}, {0, -1/3 + 1/64, -1/6 - 1/128}},
  {"near_side_outside", 0, {-2/5, 1/3, 3/4}, {0, -1/3 - 1/64, -1/6 + 1/128}},
  {"exterior_oblique", 0, {1, 0, 0}, {3/4, 1/2, 1/3}},
  {"far_reverse", 0, {1/2, -1/4, -1}, {-1, 3/4, 5/4}},
  {"shape1_interior_centroid", 1, {-1/3, 2/5, -3/4},
    altTetraCenter},
  {"shape1_near_base_inside", 1, {2/3, -1/2, 3/4},
    altTetraBaseCenter + (altTetraCenter - altTetraBaseCenter)/128},
  {"shape1_near_base_outside", 1, {2/3, -1/2, 3/4},
    altTetraBaseCenter - (altTetraCenter - altTetraBaseCenter)/128},
  {"shape1_near_edge_outside", 1, {-3/4, 1/3, 1/2},
    altTetraEdge + (altTetraBaseCenter - altTetraEdge)/32 -
      (altTetraCenter - altTetraBaseCenter)/64},
  {"shape1_far_weak", 1, {1/256, -1/384, -1/128}, {2, -3, 4}}
};

If[!TrueQ[polyInsideQ[polyTetraVertices, polyTetraFaces, {0, -1/8, -1/4}]] ||
   TrueQ[polyInsideQ[polyTetraVertices, polyTetraFaces, {0, 0, 1}]],
  Print["FAIL tetrahedron outward-face orientation"]; Exit[1]
];
If[!TrueQ[polyOutwardQ[polyTetraVerticesAlt, polyTetraFaces]],
  Print["FAIL alternate tetrahedron outward-face orientation"]; Exit[1]
];
polyGenerate["tetrahedron_magnet", tetraCases, polyTetraGeometries,
  polyClosedB[polyTetraGeometries[[#1 + 1, 1]],
    polyTetraGeometries[[#1 + 1, 2]], #2, #3, #4] &];
