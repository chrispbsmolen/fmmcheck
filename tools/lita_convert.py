"""Read a lita .npz scheme (rows of U, V, W over flattened N*N matrices, C_k = sum T_ijk A_i B_j), find the
flattening convention by a floating point screen, and write it in the catalogue tensor format for exact checking."""
import sys, json, itertools
import numpy as np
from fractions import Fraction
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from catfmt import write_full
path, out = sys.argv[1], sys.argv[2]
z = np.load(path, allow_pickle=False); meta = json.loads(z['metadata_json'].item()); N = meta['tensor'][0]; R = meta['rank']
def rows(name):
    ip, ix, nu, de = z[name + '_indptr'], z[name + '_indices'], z[name + '_numerators'], z[name + '_denominators']
    return [list(zip(ix[ip[q]:ip[q + 1]].tolist(), nu[ip[q]:ip[q + 1]].tolist(), de[ip[q]:ip[q + 1]].tolist())) for q in range(R)]
Ur, Vr, Wr = rows('u'), rows('v'), rows('w')
if '--corrupt' in sys.argv:   # same change as tools/lita_exact.py --corrupt: first U coefficient + 1
    i0, n0, d0 = Ur[0][0]; Ur[0][0] = (i0, n0 + d0, d0)
def dense(rs):
    M = np.zeros((R, N * N))
    for q, r in enumerate(rs):
        for i, n, d in r: M[q, i] += n / d
    return M
Uf, Vf, Wf = dense(Ur), dense(Vr), dense(Wr)
rng = np.random.default_rng(3); A = rng.normal(size=(N, N)); B = rng.normal(size=(N, N)); C = A @ B
best = None
for oa, ob, oc in itertools.product('rc', repeat=3):
    va = A.reshape(-1) if oa == 'r' else A.T.reshape(-1); vb = B.reshape(-1) if ob == 'r' else B.T.reshape(-1)
    ck = Wf.T @ ((Uf @ va) * (Vf @ vb)); vc = C.reshape(-1) if oc == 'r' else C.T.reshape(-1)
    err = np.abs(ck - vc).max(); print(oa, ob, oc, '%.3g' % err)
    if best is None or err < best[0]: best = (err, oa, ob, oc)
err, oa, ob, oc = best; print('convention', oa, ob, oc, 'max error %.3g' % err, 'N', N, 'rank', R, meta)
def idx(i, o): return (i // N, i % N) if o == 'r' else (i % N, i // N)
def fmt(x): return str(x.numerator) if x.denominator == 1 else '%d/%d' % (x.numerator, x.denominator)
def mat(entries, o, transpose=False):
    M = [[Fraction(0)] * N for _ in range(N)]
    for i, n, d in entries:
        a, b = idx(i, o)
        if transpose: a, b = b, a
        M[a][b] += Fraction(n, d)
    return 'Matrix(%d, %d, [%s])' % (N, N, ','.join('[' + ','.join(fmt(v) for v in r) + ']' for r in M))
# catalogue W is Z x X with A.B = sum (U.A)(V.B) W^T, so W[l][a] = coefficient of C[a][l]
trs = ['Triad([%s, %s, %s])' % (mat(Ur[q], oa), mat(Vr[q], ob), mat(Wr[q], oc, transpose=True)) for q in range(R)]
write_full(out, N, N, N, trs); print('wrote', out)
