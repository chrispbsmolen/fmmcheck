"""Rebuild a catalogue scheme from its stated construction:
    (<outer> - g) (x) <inner>  +  sum of n_k <group product>
i.e. the tensor product of an outer and an inner scheme, where some outer triads that share one factor (up to a
scalar) are grouped, and each group of size k sharing V (resp. U) together with the inner scheme is replaced by a
scheme for <k,1,1> (x) inner = <k X2, Y2, Z2> (resp. <1,1,k> (x) inner = <X2, Y2, k Z2>).
Conventions are those of the catalogue files: A.B = sum Tr(U^T A) Tr(V^T B) W^T, U X x Y, V Y x Z, W Z x X.
Correctness is NOT assumed: the output is verified exactly by verify_fmm.py."""
import sys, bz2, json, itertools
from fractions import Fraction
import os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify_fmm import parse_scheme
def load(path):
    t, _, alg = parse_scheme(path); assert not alg; return [[[[Fraction(v) for v in row] for row in M] for M in tr] for tr in t]
def T(M): return [list(r) for r in zip(*M)]
def kron(P, Q):
    return [[P[i][j] * Q[k][l] for j in range(len(P[0])) for l in range(len(Q[0]))] for i in range(len(P)) for k in range(len(Q))]
def add(P, Q): return [[a + b for a, b in zip(r, s)] for r, s in zip(P, Q)]
def scale(c, P): return [[c * a for a in r] for r in P]
def dims(s): U, V, W = s[0]; return len(U), len(U[0]), len(V[0])
def orient(s, X, Y, Z):
    """return the scheme transformed by the format symmetries to <X,Y,Z>"""
    cands = []
    for k in range(3):
        c = [tuple(t[(i + k) % 3] for i in range(3)) for t in s]       # cyclic: <X,Y,Z> -> <Y,Z,X>
        cands.append(c)
        cands.append([(T(w), T(v), T(u)) for u, v, w in c])              # transpose: Tr(ABC) = Tr(C^T B^T A^T)
    for c in cands:
        if dims(c) == (X, Y, Z): return c
    raise ValueError('no orientation gives %dx%dx%d' % (X, Y, Z))
def key(M):
    flat = [v for r in M for v in r]; nz = [v for v in flat if v]
    if not nz: return None, None
    s0 = nz[0]; return tuple(v / s0 for v in flat), s0
def build(outer, inner, groups):
    """groups: list of (k, shared, count, scheme) with shared in 'V' (group spans <k,1,1>) or 'U' (<1,1,k>)"""
    X2, Y2, Z2 = dims(inner)
    used = set(); out = []; plan = []
    for k, shared, count, gs in groups:
        idx = 0 if shared == 'U' else 1
        classes = {}
        for i, t in enumerate(outer):
            if i in used: continue
            kk, s0 = key(t[idx]); classes.setdefault(kk, []).append((i, s0))
        made = 0
        for kk, members in sorted(classes.items(), key=lambda x: -len(x[1])):
            while len(members) >= k and made < count:
                grp, members = members[:k], members[k:]
                for i, _ in grp: used.add(i)
                plan.append((k, shared, gs, grp)); made += 1
        if made < count: raise ValueError('only %d groups of %d sharing %s, need %d' % (made, k, shared, count))
    for i, (u, v, w) in enumerate(outer):
        if i in used: continue
        for x, y, z in inner: out.append((kron(u, x), kron(v, y), kron(w, z)))
    for k, shared, gs, grp in plan:
        base = outer[grp[0][0]]
        for x, y, z in gs:
            if shared == 'V':      # <k,1,1> (x) inner = <k X2, Y2, Z2>; outer triads u_p (x) v (x) w_p, v_p = s_p v
                v = base[1]
                U = None; W = None
                for p, (i, s_p) in enumerate(grp):
                    u_p = scale(s_p / grp[0][1], outer[i][0]); w_p = outer[i][2]
                    xp = x[p * X2:(p + 1) * X2]; zp = [r[p * X2:(p + 1) * X2] for r in z]
                    U = kron(u_p, xp) if U is None else add(U, kron(u_p, xp)); W = kron(w_p, zp) if W is None else add(W, kron(w_p, zp))
                out.append((U, kron(v, y), W))
            else:                  # <1,1,k> (x) inner = <X2, Y2, k Z2>; outer triads u (x) v_q (x) w_q, u_q = s_q u
                u = base[0]
                V = None; W = None
                for q, (i, s_q) in enumerate(grp):
                    v_q = scale(s_q / grp[0][1], outer[i][1]); w_q = outer[i][2]
                    yq = [r[q * Z2:(q + 1) * Z2] for r in y]; zq = z[q * Z2:(q + 1) * Z2]
                    V = kron(v_q, yq) if V is None else add(V, kron(v_q, yq)); W = kron(w_q, zq) if W is None else add(W, kron(w_q, zq))
                out.append((kron(u, x), V, W))
    return out
def fmt(x): return str(x.numerator) if x.denominator == 1 else '%d/%d' % (x.numerator, x.denominator)
def mat(M): return 'Matrix(%d, %d, [%s])' % (len(M), len(M[0]), ','.join('[' + ','.join(fmt(v) for v in r) + ']' for r in M))
def write(s, path):
    with bz2.open(path, 'wt') as f: f.write('Tensor:=TriadSet([' + ', '.join('Triad([%s, %s, %s])' % tuple(mat(M) for M in t) for t in s) + ']):\n')
