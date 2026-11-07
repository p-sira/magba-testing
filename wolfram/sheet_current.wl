(* Independent tetrahedral surface-current Biot-Savart references.
   Run: wolframscript -file testing/wolfram/sheet_current.wl *)
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "common.wl"}]];
Get[FileNameJoin[{DirectoryName[ExpandFileName[$InputFileName]], "surface_current_common.wl"}]];

ClearAll[sheetField];
sheetVertices = {
  {-1/2, -1/2, -1/2}, {1/2, -1/2, -1/2},
  {0, 1/2, -1/2}, {0, 0, 1/2}
};
sheetFaces = {{1, 3, 2}, {1, 2, 4}, {2, 3, 4}, {1, 4, 3}};
sheetField[k_, p_, wp_] := Total[
  hpSurfaceCurrentField[k, p, sheetVertices[[#]], wp] & /@ sheetFaces
];

sheetCases = {
  {"outside_top", {1, 0, 0}, {0, 0, 3/2}},
  {"outside_side", {0, 1, 0}, {3/2, 1/5, 1/6}},
  {"inside_mixed", {1/2, -3/4, 2/3}, {0, 0, -1/8}},
  {"inside_offcenter", {2/3, 1/2, -1/4}, {1/16, -1/16, -1/4}},
  {"outside_below", {0, 0, 5/4}, {0, 0, -3/2}},
  {"near_face_outside", {1, -2/3, 1/4}, {0, -1/4, -1/2 - 1/64}},
  {"near_face_inside", {-1/2, 1, -3/4}, {0, -1/4, -1/2 + 1/64}},
  {"near_top_vertex", {1, 1/3, -1/5}, {0, 0, 1/2 + 1/32}},
  {"near_bottom_edge", {3/4, 1, 0},
    {1/64, -1/2 + 1/32, -1/2 + 1/32}},
  {"far_oblique", {1/3, -2/5, 4/7}, {3, 4, 5}}
};

sheetRows = Table[
  Module[{label = item[[1]], k = item[[2]], p = item[[3]], field},
    field = hpConfirmed["sheet_current/" <> label, sheetField[k, p, #] &];
    Join[{label}, hpASCII /@ k, hpASCII /@ p, hpASCII /@ field]
  ], {item, sheetCases}
];
Print["SheetCurrent four-face tangent projections checked."];
hpWriteCSV["sheet_current",
  {"case", "kx", "ky", "kz", "x", "y", "z", "bx", "by", "bz"}, sheetRows];
