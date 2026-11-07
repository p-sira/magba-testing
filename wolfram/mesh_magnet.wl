(* Independent high-precision closed MeshMagnet surface-charge reference.
   Mesh currently consists of the four outward tetrahedron triangles.
   Run: wolframscript -file testing/wolfram/mesh_magnet.wl *)
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "common.wl"}]];
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "polyhedron_common.wl"}]];

meshVertices = polyTetraVertices;
meshFaces = polyTetraFaces;

meshCases = {
  {"interior_centroid", {0, 0, 1}, {0, -1/8, -1/4}},
  {"interior_oblique", {-1/3, 3/5, 1/2}, {0, -1/16, -1/8}},
  {"near_base_inside", {1/2, 1/3, -2/5}, {0, -1/6, -63/128}},
  {"near_base_outside", {1/2, 1/3, -2/5}, {0, -1/6, -65/128}},
  {"near_side_inside", {1/4, -3/5, 2/3}, {0, -1/3 + 1/64, -1/6 - 1/128}},
  {"near_side_outside", {1/4, -3/5, 2/3}, {0, -1/3 - 1/64, -1/6 + 1/128}},
  {"exterior_oblique", {-1/2, 0, 1}, {3/4, 1/2, 1/3}},
  {"far_reverse", {2/5, -1/2, -3/4}, {-1, 3/4, 5/4}}
};

If[!TrueQ[polyInsideQ[meshVertices, meshFaces, {0, -1/8, -1/4}]] ||
   TrueQ[polyInsideQ[meshVertices, meshFaces, {0, 0, 1}]],
  Print["FAIL mesh outward-face orientation"]; Exit[1]
];
polyGenerate["mesh_magnet", meshCases,
  polyClosedB[meshVertices, meshFaces, #1, #2, #3] &];
