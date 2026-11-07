(* Independent high-precision closed MeshMagnet surface-charge reference.
   Shapes 0 and 1 are tetrahedral meshes; shape 2 is a triangular bipyramid.
   Run: wolframscript -file testing/wolfram/mesh_magnet.wl *)
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "common.wl"}]];
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "polyhedron_common.wl"}]];

meshVertices = polyTetraVertices;
meshFaces = polyTetraFaces;
altMesh = polyTetraVerticesAlt;
altMeshCenter = Total[altMesh]/4;
altMeshSideCenter = Total[altMesh[[polyTetraFaces[[2]] + 1]]]/3;
meshBipyramid = polyBipyramidVertices;
meshBipyramidCenter = Total[meshBipyramid]/Length[meshBipyramid];
meshBipyramidTopCenter = Total[meshBipyramid[[polyBipyramidFaces[[1]] + 1]]]/3;
meshBipyramidBottomCenter =
  Total[meshBipyramid[[polyBipyramidFaces[[4]] + 1]]]/3;

meshCases = {
  {"interior_centroid", 0, {0, 0, 1}, {0, -1/8, -1/4}},
  {"interior_oblique", 0, {-1/3, 3/5, 1/2}, {0, -1/16, -1/8}},
  {"near_base_inside", 0, {1/2, 1/3, -2/5}, {0, -1/6, -63/128}},
  {"near_base_outside", 0, {1/2, 1/3, -2/5}, {0, -1/6, -65/128}},
  {"near_side_inside", 0, {1/4, -3/5, 2/3}, {0, -1/3 + 1/64, -1/6 - 1/128}},
  {"near_side_outside", 0, {1/4, -3/5, 2/3}, {0, -1/3 - 1/64, -1/6 + 1/128}},
  {"exterior_oblique", 0, {-1/2, 0, 1}, {3/4, 1/2, 1/3}},
  {"far_reverse", 0, {2/5, -1/2, -3/4}, {-1, 3/4, 5/4}},
  {"shape1_interior", 1, {-2/5, 1/3, 3/4}, altMeshCenter},
  {"shape1_near_side_inside", 1, {1/2, -2/3, 3/5},
    altMeshSideCenter + (altMeshCenter - altMeshSideCenter)/128},
  {"shape1_near_side_outside", 1, {1/2, -2/3, 3/5},
    altMeshSideCenter - (altMeshCenter - altMeshSideCenter)/128},
  {"shape1_far_weak", 1, {-1/512, 1/256, -1/384}, {-2, 3/2, 5/2}},
  {"shape2_interior", 2, {1/3, -1/2, 4/5}, meshBipyramidCenter},
  {"shape2_near_top_outside", 2, {-3/4, 2/5, -1/2},
    meshBipyramidTopCenter -
      (meshBipyramidCenter - meshBipyramidTopCenter)/128},
  {"shape2_near_bottom_inside", 2, {2/3, -1/4, 1/2},
    meshBipyramidBottomCenter +
      (meshBipyramidCenter - meshBipyramidBottomCenter)/128}
};

If[!TrueQ[polyInsideQ[meshVertices, meshFaces, {0, -1/8, -1/4}]] ||
   TrueQ[polyInsideQ[meshVertices, meshFaces, {0, 0, 1}]],
  Print["FAIL mesh outward-face orientation"]; Exit[1]
];
If[!TrueQ[polyOutwardQ[altMesh, polyTetraFaces]] ||
   !TrueQ[polyOutwardQ[meshBipyramid, polyBipyramidFaces]],
  Print["FAIL alternate mesh outward-face orientation"]; Exit[1]
];
polyGenerate["mesh_magnet", meshCases, polyMeshGeometries,
  polyClosedB[polyMeshGeometries[[#1 + 1, 1]],
    polyMeshGeometries[[#1 + 1, 2]], #2, #3, #4] &];
