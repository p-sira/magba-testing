(* Shared, Wolfram-only helpers for independent magnetic-field references.
   Exact source input values must be integers or rationals. *)

ClearAll[hpASCII, hpConfirmed, hpWriteCSV];
hpMu0 = 125663706127/10^17;
hpOutputDirectory = FileNameJoin[{
  DirectoryName[ExpandFileName[$InputFileName]], "..", "data", "high-precision"
}];

hpASCII[value_] := Module[{digits, exponent, sign},
  If[!NumericQ[value],
    Print["FAIL nonnumeric output: ", InputForm[value]]; Exit[1]
  ];
  (* Round away symmetry residuals many orders beneath f64 resolution. *)
  If[PossibleZeroQ[value] || TrueQ[Abs[value] < 10^-40], Return["0"]];
  {digits, exponent} = RealDigits[Abs[value], 10, 26];
  If[!VectorQ[digits, IntegerQ],
    Print["FAIL nonnumeric decimal digits: ", InputForm[value]]; Exit[1]
  ];
  sign = If[TrueQ[value < 0], "-", ""];
  sign <> ToString[First[digits]] <> "." <>
    StringJoin[ToString /@ Rest[digits]] <> "e" <> ToString[exponent - 1]
];

hpConfirmed[label_String, evaluator_] := Module[{low, high, error, scale},
  low = Check[evaluator[40], $Failed, {
    NIntegrate::eincr, NIntegrate::ncvb, NIntegrate::slwcon,
    NIntegrate::inumr, NIntegrate::precw
  }];
  high = Check[evaluator[60], $Failed, {
    NIntegrate::eincr, NIntegrate::ncvb, NIntegrate::slwcon,
    NIntegrate::inumr, NIntegrate::precw
  }];
  If[!VectorQ[low, NumericQ] || !VectorQ[high, NumericQ] ||
     Length[low] != 3 || Length[high] != 3,
    Print["FAIL ", label, ": nonnumeric or nonconverged field"]; Exit[1]
  ];
  error = Max[Abs[high - low]];
  scale = Max[10^-18, Max[Abs[high]]];
  If[!TrueQ[error < 10^-20 scale],
    Print["FAIL ", label, ": 40/60-digit disagreement ",
      hpASCII[error/scale]]; Exit[1]
  ];
  Print[label, ": 40/60-digit agreement < 1e-20 scaled"];
  high
];

hpWriteCSV[name_String, headers_List, rows_List] := Module[{path},
  If[!DirectoryQ[hpOutputDirectory],
    CreateDirectory[hpOutputDirectory, CreateIntermediateDirectories -> True]
  ];
  path = FileNameJoin[{hpOutputDirectory, name <> ".csv"}];
  Export[path, Prepend[rows, headers], "CSV"];
  Print["Wrote ", Length[rows], " cases to ", path];
];
