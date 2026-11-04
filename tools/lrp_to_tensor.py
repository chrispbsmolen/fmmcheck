"""Convert the catalogue's LRP representation (L: r x XY, R: r x YZ, P: XZ x r; products
m_k = (L_k . vec A)(R_k . vec B), vec(C) = P m) into the catalogue's tensor format
(TriadSet of Triad([U (X x Y), V (Y x Z), W (Z x X)]) with A.B = sum Tr(U^T A) Tr(V^T B) W^T).
The vec orderings are found by lrp_probe.py (floating point screen over the 8 combinations) and
the converted file is then verified exactly by verify_fmm.py; nothing is assumed."""
import re, sys, bz2
from fractions import Fraction
X, Y, Z = map(int, sys.argv[1:4]); src = sys.argv[4]; out = sys.argv[5]; oa, ob, oc = sys.argv[6]  # vec order of A, B, C: r (row-major) or c (column-major)
s = bz2.open(src, 'rt').read()
mats = []
for m in re.finditer(r'Matrix\((\d+), (\d+), \[\[(.*?)\]\]\)', s, re.S):
    r, c = int(m.group(1)), int(m.group(2))
    rows = m.group(3).replace('\n', '').split('],[')
    M = [[Fraction(x) for x in row.split(',')] for row in rows]
    assert len(M) == r and all(len(x) == c for x in M); mats.append(M)
L, R, P = mats
rank = len(L); assert len(R) == rank and len(P[0]) == rank and len(L[0]) == X * Y and len(R[0]) == Y * Z and len(P) == X * Z
def unvec(v, n, m, o):  # vector of length n*m -> n x m matrix
    return [[v[i * m + j] if o == 'r' else v[j * n + i] for j in range(m)] for i in range(n)]
def fmt(x): return str(x.numerator) if x.denominator == 1 else '%d/%d' % (x.numerator, x.denominator)
def mat(M): return 'Matrix(%d, %d, [%s])' % (len(M), len(M[0]), ','.join('[' + ','.join(fmt(x) for x in row) + ']' for row in M))
triads = []
for k in range(rank):
    U = unvec(L[k], X, Y, oa); V = unvec(R[k], Y, Z, ob)
    Cc = unvec([P[i][k] for i in range(X * Z)], X, Z, oc)          # contribution of m_k to C (X x Z)
    W = [[Cc[i][j] for i in range(X)] for j in range(Z)]          # W = Cc^T, Z x X
    triads.append('Triad([%s, %s, %s])' % (mat(U), mat(V), mat(W)))
with bz2.open(out, 'wt') as f:
    f.write('Tensor:=TriadSet([' + ', '.join(triads) + ']):\n')
print('wrote', out, 'rank', rank)
