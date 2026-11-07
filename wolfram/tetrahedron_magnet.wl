(* Independent high-precision closed TetrahedronMagnet surface-charge reference.
   Run: wolframscript -file testing/wolfram/tetrahedron_magnet.wl *)
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "common.wl"}]];
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "polyhedron_common.wl"}]];

tetraCases = {
  {"interior_centroid", {0, 0, 1}, {0, -1/8, -1/4}},
  {"interior_vector", {1/3, -2/5, 3/4}, {1/16, -1/8, -1/4}},
  {"near_base_inside", {0, 0, 1}, {0, -1/6, -63/128}},
  {"near_base_outside", {0, 0, 1}, {0, -1/6, -65/128}},
  {"near_side_inside", {-2/5, 1/3, 3/4}, {0, -1/3 + 1/64, -1/6 - 1/128}},
  {"near_side_outside", {-2/5, 1/3, 3/4}, {0, -1/3 - 1/64, -1/6 + 1/128}},
  {"exterior_oblique", {1, 0, 0}, {3/4, 1/2, 1/3}},
  {"far_reverse", {1/2, -1/4, -1}, {-1, 3/4, 5/4}}
};

If[!TrueQ[polyInsideQ[polyTetraVertices, polyTetraFaces, {0, -1/8, -1/4}]] ||
   TrueQ[polyInsideQ[polyTetraVertices, polyTetraFaces, {0, 0, 1}]],
  Print["FAIL tetrahedron outward-face orientation"]; Exit[1]
];
polyGenerate["tetrahedron_magnet", tetraCases,
  polyClosedB[polyTetraVertices, polyTetraFaces, #1, #2, #3] &];
