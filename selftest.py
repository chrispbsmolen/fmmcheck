"""Tests for verify_fmm.py.  Run from the repository root:  python3 selftest.py

Part A needs nothing but this repository.  It builds schemes in a temporary folder (the naive algorithm,
Strassen, tensor products, the same schemes after random rational basis changes, and parametric entries),
corrupts some of them in the ways a checker could miss (a coefficient off by 1, off by exactly a CRT prime
or a product of two, a dropped triad, an error in the last triad), and compares verify() under every
arithmetic path (int64, CRT, dense, dense with several block sizes) with a brute-force oracle written here
that evaluates every Brent coefficient with Python Fractions.  It also checks Maple line wrapping inside a
token and backslash continuations.

Part B runs only if data/tensors/ holds catalogue files (python3 fetch_catalogue.py): for a fixed random
sample of small formats it compares verify() with the oracle, and checks a corrupted copy of each.

Exit status 0 only if every test passes."""
import os, sys, bz2, json, random, itertools, subprocess, tempfile
from fractions import Fraction as F
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from verify_fmm import parse_scheme
TMP = tempfile.mkdtemp(prefix='fmmcheck_selftest_'); P1 = 1073741789; P2 = 1073741783
failures = 0
def E(r, c, i, j): M = [[0] * c for _ in range(r)]; M[i][j] = 1; return M
def naive(X, Y, Z): return [(E(X, Y, a, b), E(Y, Z, b, d), E(Z, X, d, a)) for a in range(X) for b in range(Y) for d in range(Z)]
def strassen():
    m = lambda l: [[l[0], l[1]], [l[2], l[3]]]
    pr = [((1,0,0,1),(1,0,0,1),(1,0,0,1)), ((0,0,1,1),(1,0,0,0),(0,0,1,-1)), ((1,0,0,0),(0,1,0,-1),(0,1,0,1)), ((0,0,0,1),(-1,0,1,0),(1,0,1,0)),
          ((1,1,0,0),(0,0,0,1),(-1,1,0,0)), ((-1,0,1,0),(1,1,0,0),(0,0,0,1)), ((0,1,0,-1),(0,0,1,1),(1,0,0,0))]
    return [(m(u), m(v), [list(r) for r in zip(*m(c))]) for u, v, c in pr]
def kron(A, B): return [[A[i // len(B)][j // len(B[0])] * B[i % len(B)][j % len(B[0])] for j in range(len(A[0]) * len(B[0]))] for i in range(len(A) * len(B))]
def tprod(S1, S2): return [tuple(kron(a, b) for a, b in zip(t1, t2)) for t1 in S1 for t2 in S2]
def mm(A, B): return [[sum(F(A[i][k]) * B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def inv(M):
    n = len(M); A = [[F(x) for x in r] + [F(int(i == j)) for j in range(n)] for i, r in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c]); A[c], A[p] = A[p], A[c]; A[c] = [x / A[c][c] for x in A[c]]
        for r in range(n):
            if r != c and A[r][c]: A[r] = [x - A[r][c] * y for x, y in zip(A[r], A[c])]
    return [r[n:] for r in A]
def T(M): return [list(r) for r in zip(*M)]
def basis_change(S, X, Y, Z, rng):
    """A -> P A Q^-1, B -> Q B R^-1, C -> P C R^-1: U' = P^T U Q^-T, V' = Q^T V R^-T, W' = R^T... (W = C-coefficients^T)"""
    def rnd(n):
        while True:
            M = [[F(rng.randint(-3, 3), rng.choice([1, 1, 2, 3])) for _ in range(n)] for _ in range(n)]
            try: inv(M); return M
            except StopIteration: pass
    Pm, Qm, Rm = rnd(X), rnd(Y), rnd(Z); Pi, Qi, Ri = inv(Pm), inv(Qm), inv(Rm)
    out = []
    for U, V, W in S:
        N = T(W)                                           # coefficient matrix of C (X x Z)
        out.append((mm(mm(T(Pm), U), T(Qi)), mm(mm(T(Qm), V), T(Ri)), T(mm(mm(Pi, N), Rm))))
    return out
def oracle(S, X, Y, Z):
    """brute force: every Brent coefficient with Fractions; returns (wrong, missing) with the verifier's meaning"""
    acc = {}; 
    for U, V, W in S:
        for a, b, c, d, f, e in itertools.product(range(X), range(Y), range(Y), range(Z), range(Z), range(X)):
            u, v, w = U[a][b], V[c][d], W[f][e]
            if u and v and w: acc[(a, b, c, d, e, f)] = acc.get((a, b, c, d, e, f), 0) + F(u) * v * w
    wrong = sum(1 for (a, b, c, d, e, f), s in acc.items() if s != (1 if (b == c and a == e and d == f) else 0))
    missing = sum(1 for a in range(X) for b in range(Y) for d in range(Z) if (a, b, b, d, a, d) not in acc)
    return wrong, missing
def fmt(v): v = F(v); return str(v.numerator) if v.denominator == 1 else '%d/%d' % (v.numerator, v.denominator)
def write(S, name, wrap=None):
    mat = lambda M: 'Matrix(%d, %d, [%s])' % (len(M), len(M[0]), ','.join('[' + ','.join(fmt(x) for x in r) + ']' for r in M))
    body = 'Tensor:=TriadSet([' + ', '.join('Triad([%s, %s, %s])' % tuple(mat(M) for M in t) for t in S) + ']):\n'
    if wrap: body = wrap(body)
    p = os.path.join(TMP, name + '.mpl.bz2'); bz2.open(p, 'wt').write(body); return p
def run(p, X, Y, Z, r, env):
    e = {k: v for k, v in os.environ.items() if not k.startswith('FMM_')}; e.update(env)
    o = subprocess.run([sys.executable, '-c', 'import sys,json;sys.path.insert(0,%r);from verify_fmm import verify;print(json.dumps(verify(sys.argv[1],*map(int,sys.argv[2:]))))' % HERE, p, str(X), str(Y), str(Z), str(r)], env=e, capture_output=True, text=True)
    return json.loads(o.stdout.strip().splitlines()[-1]) if o.stdout.strip() else dict(verdict='CRASH', note=o.stderr[-200:])
MODES = [('default', {}), ('crt', {'FMM_FORCE_CRT': '1'}), ('dense', {'FMM_FORCE_DENSE': '1'}), ('dense_crt', {'FMM_FORCE_DENSE': '1', 'FMM_FORCE_CRT': '1'}),
         ('dense_blk1', {'FMM_FORCE_DENSE': '1', 'FMM_DENSE_BLOCK': '1'}), ('dense_crt_blk3', {'FMM_FORCE_DENSE': '1', 'FMM_FORCE_CRT': '1', 'FMM_DENSE_BLOCK': '3'})]
def check(label, S, X, Y, Z, stated=None, wrap=None, must_contain=None):
    global failures
    stated = len(S) if stated is None else stated
    w, m = oracle(S, X, Y, Z); expect = 'FAIL' if (w or m) else ('PASS' if len(S) == stated else 'PASS_RANK_MISMATCH')
    p = write(S, label, wrap); bad = []
    if must_contain and must_contain not in bz2.open(p, 'rt').read(): bad.append('wrapping not applied')
    for mode, env in MODES:
        r = run(p, X, Y, Z, stated, env)
        got = (r['verdict'], r.get('wrong_entries') or 0, r.get('missing_entries') or 0)
        if got != (expect, w, m) or r.get('rank') != len(S): bad.append('%s:%s' % (mode, got))
    print('%-34s expected %-18s wrong %-4d missing %-3d %s' % (label, expect, w, m, 'ok' if not bad else 'PROBLEM ' + ' '.join(bad)), flush=True)
    if bad: failures += 1
def corrupt(S, k, i, j, delta, which=0):
    S2 = [tuple([list(map(list, M)) for M in t]) for t in S]; M = S2[k][which]
    M[i][j] = F(M[i][j]) + delta; return S2
print('Part A: constructed schemes')
rng = random.Random(20261104)
st = strassen(); st2 = tprod(st, st); stn = tprod(st, naive(1, 2, 2))
check('naive_2x3x4', naive(2, 3, 4), 2, 3, 4)
check('strassen', st, 2, 2, 2)
check('strassen_x_strassen', st2, 4, 4, 4)
check('strassen_x_naive_1x2x2', stn, 2, 4, 4)
bst = basis_change(st, 2, 2, 2, rng); check('strassen_basis_change', bst, 2, 2, 2)
bst2 = basis_change(st2, 4, 4, 4, rng); check('strassen2_basis_change', bst2, 4, 4, 4)
check('rank_mismatch', st, 2, 2, 2, stated=8)
check('plus1_first', corrupt(st, 0, 0, 0, 1), 2, 2, 2)
check('plus1_last_triad', corrupt(st2, len(st2) - 1, 3, 3, 1, 2), 4, 4, 4)
check('plus_P1', corrupt(st, 2, 0, 0, P1), 2, 2, 2)
check('plus_P1P2', corrupt(st, 4, 0, 1, P1 * P2, 1), 2, 2, 2)
check('plus_1_over_P1', corrupt(bst, 3, 1, 1, F(1, P1)), 2, 2, 2)
check('dropped_triad', st2[:-1], 4, 4, 4)
check('empty_triad_added', st + [([[0, 0], [0, 0]], st[0][1], st[0][2])], 2, 2, 2, stated=7)
check('huge_denominator', [(U, [[F(x) * F(1, 2 ** 70) for x in r] for r in V], [[F(x) * 2 ** 70 for x in r] for r in W]) for U, V, W in st], 2, 2, 2)
# Maple line wrapping inside a token and backslash continuation
def wrap_inside(b):   # split a numerator, a denominator and the word Matrix across lines
    i = b.index('/'); b = b[:i - 1] + '\n' + b[i - 1:]; j = b.index('/', i + 3); b = b[:j + 1] + '\n' + b[j + 1:]
    return b.replace('Matrix(', 'Mat\nrix(', 2)
check('newline_inside_tokens', bst, 2, 2, 2, wrap=wrap_inside, must_contain='Mat\nrix(')
check('backslash_continuation', st2, 4, 4, 4, wrap=lambda b: b.replace(',1,', ',1\\\n,', 5), must_contain='\\\n')
# parametric entries: Strassen with M1 scaled by l in U and 1/l in W
S_par = write(st, 'par_base')
txt = bz2.open(S_par, 'rt').read().replace('Triad([Matrix(2, 2, [[1,0],[0,1]])', 'Triad([Matrix(2, 2, [[l,0],[0,l]])', 1)
i = txt.find('Matrix(2, 2, [[1,0],[0,1]])])'); txt = txt[:i] + 'Matrix(2, 2, [[1/l,0],[0,1/l]])])' + txt[i + len('Matrix(2, 2, [[1,0],[0,1]])])'):]
pp = os.path.join(TMP, 'parametric.mpl.bz2'); bz2.open(pp, 'wt').write(txt)
r = run(pp, 2, 2, 2, 7, {}); ok = r['verdict'] == 'PASS' and r.get('ring') == 'parametric'
print('%-34s expected %-18s %s' % ('parametric_l', 'PASS (parametric)', 'ok' if ok else 'PROBLEM %s' % r)); failures += not ok
print('\nPart B: catalogue files')
D = os.path.join(HERE, 'data')
if not os.path.isdir(D + '/tensors') or not os.listdir(D + '/tensors'):
    print('skipped: data/tensors/ is empty (run python3 fetch_catalogue.py first)')
else:
    Ent = json.load(open(D + '/entries.json'))
    small = sorted(f for f, e in Ent.items() if e['X'] * e['Y'] * e['Z'] <= 60 and os.path.exists(D + '/tensors/%s_tensor.mpl.bz2' % f))
    random.seed(7)
    for f in random.sample(small, min(15, len(small))):
        e = Ent[f]; S, _, alg = parse_scheme(D + '/tensors/%s_tensor.mpl.bz2' % f)
        if alg: print('%-34s skipped (parametric)' % f); continue
        check('catalogue_' + f, S, e['X'], e['Y'], e['Z'], stated=e['rank'])
        k = next(i for i, t in enumerate(S) if any(v for row in t[0] for v in row))
        a, b = next((i, j) for i, row in enumerate(S[k][0]) for j, v in enumerate(row) if v)
        check('catalogue_' + f + '_corrupted', corrupt(S, k, a, b, 1), e['X'], e['Y'], e['Z'], stated=e['rank'])
print('\n%s' % ('ALL TESTS PASS' if failures == 0 else '%d PROBLEMS' % failures))
sys.exit(1 if failures else 0)
