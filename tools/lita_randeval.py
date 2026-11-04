"""Second check of a lita .npz scheme, independent of tools/lita_exact.py: multiply random integer N x N matrices
with the scheme's formula in exact rational arithmetic and count wrong entries of A.B.
    python3 tools/lita_randeval.py SCHEME.npz"""
import sys, json, random
import numpy as np
from fractions import Fraction
z = np.load(sys.argv[1], allow_pickle=False); meta = json.loads(z['metadata_json'].item()); N = meta['tensor'][0]; R = meta['rank']
rows = lambda n: (z[n + '_indptr'], z[n + '_indices'].tolist(), z[n + '_numerators'].tolist(), z[n + '_denominators'].tolist())
(ui, ux, un, ud), (vi, vx, vn, vd), (wi, wx, wn, wd) = rows('u'), rows('v'), rows('w')
rng = random.Random(2026); A = [rng.randint(-9, 9) for _ in range(N * N)]; B = [rng.randint(-9, 9) for _ in range(N * N)]
C = [Fraction(0)] * (N * N)
for q in range(R):
    s = sum(Fraction(un[t], ud[t]) * A[ux[t]] for t in range(ui[q], ui[q + 1]))
    r = sum(Fraction(vn[t], vd[t]) * B[vx[t]] for t in range(vi[q], vi[q + 1]))
    if s and r:
        sr = s * r
        for t in range(wi[q], wi[q + 1]): C[wx[t]] += sr * Fraction(wn[t], wd[t])
true = [sum(A[a * N + b] * B[b * N + d] for b in range(N)) for a in range(N) for d in range(N)]   # row-major C = A.B
print(json.dumps(dict(scheme=sys.argv[1].split('/')[-1], entries=N * N, wrong=sum(1 for x, y in zip(C, true) if x != y))))
