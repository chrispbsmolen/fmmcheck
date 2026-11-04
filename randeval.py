"""Second, independent check: evaluate a file's own identity A.B = sum_i Tr(U_i^T A) Tr(V_i^T B) W_i^T on
random integer matrices A (X x Y) and B (Y x Z) in exact rational arithmetic, and count the wrong entries
of A.B.  Shares only the file parser with verify_fmm.py, none of the accumulation code.  A wrong entry
proves the file wrong; agreement on random inputs is evidence, not proof (verify_fmm.py is the proof).

  python3 randeval.py FILE X Y Z [trials]"""
import sys, json, random
from fractions import Fraction
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_fmm import parse_scheme
path = sys.argv[1]; X, Y, Z = map(int, sys.argv[2:5]); trials = int(sys.argv[5]) if len(sys.argv) > 5 else 1; f = os.path.basename(path)
triads, _, alg = parse_scheme(path)
assert not alg
rng = random.Random(2026)
for t in range(trials):
    A = [[rng.randint(-9, 9) for _ in range(Y)] for _ in range(X)]
    B = [[rng.randint(-9, 9) for _ in range(Z)] for _ in range(Y)]
    C = [[Fraction(0)] * X for _ in range(Z)]   # accumulates (A.B)^T, Z x X
    for U, V, W in triads:
        s = sum(U[a][b] * A[a][b] for a in range(X) for b in range(Y) if U[a][b])
        if not s: continue
        r = sum(V[c][d] * B[c][d] for c in range(Y) for d in range(Z) if V[c][d])
        if not r: continue
        sr = s * r
        for i in range(Z):
            Wi = W[i]
            for j in range(X):
                if Wi[j]: C[i][j] += sr * Wi[j]
    AB = [[sum(A[a][b] * B[b][d] for b in range(Y)) for d in range(Z)] for a in range(X)]
    bad = [(a, d) for a in range(X) for d in range(Z) if C[d][a] != AB[a][d]]
    ex = None if not bad else dict(entry=list(bad[0]), scheme_gives=str(C[bad[0][1]][bad[0][0]]), true_value=AB[bad[0][0]][bad[0][1]])
    print(json.dumps(dict(file=f, trial=t, entries=X * Z, wrong=len(bad), first_wrong=ex)), flush=True)
