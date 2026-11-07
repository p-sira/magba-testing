(* Independent tetrahedral surface-current Biot-Savart references.
   Run: wolframscript -file testing/wolfram/sheet_current.wl *)
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "common.wl"}]];
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "surface_current_common.wl"}]];

ClearAll[sheetField];
sheetVertices = {
  {
    {-1/2, -1/2, -1/2}, {1/2, -1/2, -1/2},
    {0, 1/2, -1/2}, {0, 0, 1/2}
  },
  {
    {-3/4, -1/2, -1/3}, {1/2, -1/4, -1/2},
    {1/4, 3/4, -1/4}, {-1/8, 1/8, 4/5}
  }
};
sheetFaces = {{1, 3, 2}, {1, 2, 4}, {2, 3, 4}, {1, 4, 3}};
sheetField[k_, p_, vertices_, wp_] := Total[
  hpSurfaceCurrentField[k, p, vertices[[#]], wp] & /@ sheetFaces
];
sheetBase1 = sheetVertices[[2, {1, 3, 2}]];
sheetBaseCentroid1 = Total[sheetBase1]/3;
sheetBaseNormal1 = Cross[
  sheetBase1[[2]] - sheetBase1[[1]],
  sheetBase1[[3]] - sheetBase1[[1]]
];
sheetEdgeMidpoint1 = (sheetVertices[[2, 1]] + sheetVertices[[2, 2]])/2;

sheetCases = {
  {"outside_top", 0, {1, 0, 0}, {0, 0, 3/2}},
  {"outside_side", 0, {0, 1, 0}, {3/2, 1/5, 1/6}},
  {"inside_mixed", 0, {1/2, -3/4, 2/3}, {0, 0, -1/8}},
  {"inside_offcenter", 0, {2/3, 1/2, -1/4}, {1/16, -1/16, -1/4}},
  {"outside_below", 0, {0, 0, 5/4}, {0, 0, -3/2}},
  {"near_face_outside", 0, {1, -2/3, 1/4}, {0, -1/4, -1/2 - 1/64}},
  {"near_face_inside", 0, {-1/2, 1, -3/4}, {0, -1/4, -1/2 + 1/64}},
  {"near_top_vertex", 0, {1, 1/3, -1/5}, {0, 0, 1/2 + 1/32}},
  {"near_bottom_edge", 0, {3/4, 1, 0},
    {1/64, -1/2 + 1/32, -1/2 + 1/32}},
  {"far_oblique", 0, {1/3, -2/5, 4/7}, {3, 4, 5}},
  {"shape1_near_base_outside", 1, {3/4, -2/3, 1/4},
    sheetBaseCentroid1 + sheetBaseNormal1/64},
  {"shape1_near_base_inside", 1, {-1/2, 3/4, 1/3},
    sheetBaseCentroid1 - sheetBaseNormal1/64},
  {"shape1_near_edge_inside", 1, {1, 1/2, -2/5},
    sheetEdgeMidpoint1 +
      (sheetVertices[[2, 3]] - sheetEdgeMidpoint1)/32 +
      (sheetVertices[[2, 4]] - sheetEdgeMidpoint1)/32},
  {"shape1_far_scale", 1, {-5/4, 2/3, -1/2}, {5/2, -7/3, 9/4}}
};

sheetRows = Table[
  Module[{label = item[[1]], shape = item[[2]], k = item[[3]],
          p = item[[4]], field},
    field = hpConfirmed["sheet_current/" <> label,
      sheetField[k, p, sheetVertices[[shape + 1]], #] &];
    Join[{label, ToString[shape]}, hpASCII /@ k,
      hpASCII /@ p, hpASCII /@ field]
  ], {item, sheetCases}
];
Print["SheetCurrent four-face tangent projections checked."];
hpWriteCSV["sheet_current",
  {"case", "shape", "kx", "ky", "kz", "x", "y", "z", "bx", "by", "bz"}, sheetRows];
