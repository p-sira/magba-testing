(* Independent Biot-Savart area quadrature shared by triangle and mesh fixtures.
   A current-density vector is projected onto each triangle's tangent plane. *)
ClearAll[hpSurfaceLineIntegral, hpSurfaceCurrentField];

(* On a barycentric triangle r(u,v)=a+u d+v f, 0<=u<=1-v,
   integrate (p-r)/|p-r|^3 analytically in u. *)
hpSurfaceLineIntegral[p_, a_, d_, f_, v_] := Module[
  {len = Sqrt[d.d], e, h, t0, t1, q, q2, r0, r1},
  e = d/len;
  h = p - a - v f;
  t0 = h.e;
  t1 = t0 - (1 - v) len;
  q = h - t0 e;
  q2 = q.q;
  r0 = Sqrt[q2 + t0^2];
  r1 = Sqrt[q2 + t1^2];
  (q (t0/r0 - t1/r1)/q2 + e (1/r1 - 1/r0))/len
];

hpSurfaceCurrentField[k_, p_, vertices_, wp_] := Module[
  {a = vertices[[1]], d, f, normal, tangent, jac},
  d = vertices[[2]] - a;
  f = vertices[[3]] - a;
  jac = Sqrt[Cross[d, f].Cross[d, f]];
  normal = Cross[d, f]/jac;
  tangent = k - (k.normal) normal;
  If[tangent == {0, 0, 0}, Return[{0, 0, 0}]];
  hpMu0 jac/(4 Pi) NIntegrate[
    Evaluate[Cross[tangent, hpSurfaceLineIntegral[p, a, d, f, v]]],
    {v, 0, 1}, WorkingPrecision -> wp,
    AccuracyGoal -> wp - 12, PrecisionGoal -> wp - 12,
    Method -> {"GlobalAdaptive", "SymbolicProcessing" -> 0},
    MaxRecursion -> 24
  ]
];
