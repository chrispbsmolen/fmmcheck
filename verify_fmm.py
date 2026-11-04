"""Exact verifier for the schemes in Sedoglavic's catalogue of fast matrix multiplication algorithms
(https://fmm.univ-lille.fr).

A scheme for <X x Y x Z : r> is a list of r triads (U_i, V_i, W_i), U_i an X x Y matrix, V_i a Y x Z
matrix and W_i a Z x X matrix, claimed to satisfy, for all A and B,

    A B = sum_i Tr(U_i^T A) Tr(V_i^T B) W_i^T

(the identity written at the end of every catalogue file).  Equivalently the Brent equations: for all
a, b, c, d, e, f,

    sum_i U_i[a,b] V_i[c,d] W_i[f,e] = [b = c] [a = e] [d = f].

verify(path, X, Y, Z, stated_rank) parses the file, checks the shapes, and evaluates the left-hand side
exactly.  Rational entries are cleared by the lcm of the denominators of each factor, so all arithmetic
is on integers.  Parametric entries (polynomials in free parameters) are instantiated at fixed
pseudo-random rationals (seed 12345) and the instantiated scheme is verified exactly.  The method used
is recorded in the result's 'path' field:

  numpy_int64   sparse accumulation with numpy, when a rigorous bound (sum over triads of the products of
                the l1 norms of the three factors, plus the target) is below 2^62, so nothing can overflow
  numpy_crt     the same, repeated modulo primes below 2^30 until their product exceeds the bound; equal
                residues modulo every prime then force equality (Chinese remainder theorem), and a
                difference modulo any one prime proves a wrong coefficient
  dense_int64,  for schemes with more than 2e8 terms: the coefficient tensor is held densely in blocks,
  dense_crt     with the same int64 / CRT arithmetic
  python_exact  Python integers (or sympy for non-rational algebraic entries), used when numpy is not
                applicable

Verdicts: PASS (equations hold, triad count equals the stated rank), PASS_RANK_MISMATCH (equations hold,
count differs), FAIL (wrong_entries coefficients differ among those that receive terms, missing_entries
target coefficients receive none), BADSHAPE, UNPARSED (no TriadSet in the file, or unreadable).

Environment switches used only by the tests: FMM_FORCE_CRT, FMM_FORCE_DENSE, FMM_DENSE_BLOCK.
"""
import os, re, bz2, sys, json
from fractions import Fraction
from collections import defaultdict

TOK = re.compile(r'Matrix\(\s*(\d+)\s*,\s*(\d+)\s*,\s*\[(.*?)\]\s*\)', re.S)

class Sym(str):
    """A symbolic (parametric or algebraic) entry kept as its source text."""
    pass

import ast, operator
_BIN = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul, ast.Div: operator.truediv, ast.Pow: operator.pow}

def free_symbols(s):
    return set(re.findall(r'[A-Za-z_][A-Za-z_0-9]*', s))

def eval_fraction(s, env):
    """Evaluate a Maple expression built from integers, +, -, *, /, ^ (or **) and the names in env
    exactly over Fractions.  Anything else (functions, unknown names) raises ValueError."""
    import sys
    if sys.getrecursionlimit() < 100000: sys.setrecursionlimit(100000)
    tree = ast.parse(s.replace('^', '**'), mode='eval')
    def ev(n):
        # flatten long left-associative chains of + - * so recursion depth stays small
        if isinstance(n, ast.BinOp) and type(n.op) in (ast.Add, ast.Sub, ast.Mult):
            chain = []; cur = n
            while isinstance(cur, ast.BinOp) and type(cur.op) in (ast.Add, ast.Sub, ast.Mult):
                chain.append((type(cur.op), cur.right)); cur = cur.left
            val = ev(cur)
            for op, rhs in reversed(chain): val = _BIN[op](val, ev(rhs))
            return val
        if isinstance(n, ast.Expression): return ev(n.body)
        if isinstance(n, ast.Constant) and isinstance(n.value, int) and not isinstance(n.value, bool): return Fraction(n.value)
        if isinstance(n, ast.Name):
            if n.id in env: return env[n.id]
            raise ValueError('unknown name ' + n.id)
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.USub): return -ev(n.operand)
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, ast.UAdd): return ev(n.operand)
        if isinstance(n, ast.BinOp) and type(n.op) in _BIN:
            if isinstance(n.op, ast.Pow):
                e = ev(n.right)
                if e.denominator != 1: raise ValueError('non-integer exponent')
                return ev(n.left) ** int(e)
            return _BIN[type(n.op)](ev(n.left), ev(n.right))
        raise ValueError('unsupported syntax ' + type(n).__name__)
    return ev(tree)

def parse_entry(s, algebraic):
    s = s.strip()
    try:
        return int(s)          # plain integers (the vast majority): ints behave exactly like Fractions here
    except ValueError:
        pass
    try:
        return Fraction(s)
    except Exception:
        pass
    m = re.fullmatch(r'(-?\d+)\s*/\s*(\d+)', s)
    if m: return Fraction(int(m.group(1)), int(m.group(2)))
    algebraic.append(s)
    return Sym(s)

def parse_scheme(path):
    txt = bz2.open(path, 'rt').read()
    m = re.search(r'Tensor\s*:=\s*TriadSet\(\[(.*)\]\)\s*:', txt, re.S)
    if not m: raise ValueError('no TriadSet')
    body = re.sub(r'\\\r?\n|\s+', '', m.group(1))   # Maple wraps long output lines every 1e6 characters, sometimes inside a token
    mats = []; algebraic = []
    for r, c, inner in TOK.findall(body):
        r, c = int(r), int(c)
        rows = re.findall(r'\[([^\[\]]*)\]', inner)
        M = [[parse_entry(x, algebraic) for x in row.split(',')] for row in rows]
        if len(M) != r or any(len(row) != c for row in M): raise ValueError('matrix shape mismatch %dx%d' % (r, c))
        mats.append(M)
    if len(mats) % 3: raise ValueError('matrix count not multiple of 3: %d' % len(mats))
    triads = [mats[3 * i: 3 * i + 3] for i in range(len(mats) // 3)]
    return triads, txt, algebraic

def nz(M):
    return [(i, j, v) for i, row in enumerate(M) for j, v in enumerate(row) if v != 0]

def verify(path, X, Y, Z, stated_rank):
    res = dict(X=X, Y=Y, Z=Z, stated_rank=stated_rank)
    try:
        triads, txt, algebraic = parse_scheme(path)
    except Exception as ex:
        res['verdict'] = 'UNPARSED'; res['note'] = repr(ex)[:200]; return res
    res['rank'] = len(triads)
    res['ring'] = 'algebraic' if algebraic else ('integer' if all(v.denominator == 1 for t in triads for M in t for row in M for v in row) else 'rational')
    if algebraic:
        # Parametric or algebraic entries: instantiate free symbols at random rational values
        # (exact Fractions) and verify the instantiated scheme exactly.  This proves that the
        # scheme is a valid algorithm of its rank for those parameter values, which is what
        # the catalogue's rank claim needs.  Values are drawn to avoid zero denominators.
        import random
        syms = set()
        for x in set(algebraic): syms |= free_symbols(x)
        rng = random.Random(12345); inst = None; err = None
        for attempt in range(20):
            env = {sy: Fraction(rng.randint(2, 97), rng.randint(2, 97)) for sy in sorted(syms)}   # sorted: independent of PYTHONHASHSEED
            cache = {}
            try:
                def conv(v):
                    if not isinstance(v, Sym): return v
                    if v not in cache: cache[v] = eval_fraction(v, env)
                    return cache[v]
                inst = [[[[conv(v) for v in row] for row in M] for M in t] for t in triads]
                break
            except ZeroDivisionError:
                inst = None
            except ValueError as ex:
                inst = None; err = str(ex); break
        if inst is None:
            # not a rational-parametric scheme (e.g. involves sqrt or RootOf): use sympy, slowly
            try:
                import sympy
                vals = None
                def conv2(v):
                    return sympy.nsimplify(sympy.sympify(str(v))) if isinstance(v, Sym) else sympy.Rational(v.numerator, v.denominator)
                triads = [[[[conv2(v) for v in row] for row in M] for M in t] for t in triads]
                res['ring'] = 'algebraic'; res['note'] = 'symbolic path: ' + str(err)
                algebraic_sym = True
            except Exception as ex:
                res['verdict'] = 'UNPARSED'; res['note'] = 'could not instantiate parameters: %s' % ex; return res
        else:
            triads = inst; res['ring'] = 'parametric' if syms else 'rational'; res['instantiated_symbols'] = sorted(syms)
            algebraic = []; algebraic_sym = False
    else:
        algebraic_sym = False
    # shapes
    shapes = {(len(U), len(U[0]), len(V), len(V[0]), len(W), len(W[0])) for U, V, W in triads}
    if shapes != {(X, Y, Y, Z, Z, X)}:
        res['verdict'] = 'BADSHAPE'; res['note'] = 'shapes %s expected (%d,%d,%d,%d,%d,%d)' % (sorted(shapes), X, Y, Y, Z, Z, X); return res
    # Exact integer accumulation: scale every U_i by L_U (the lcm of all denominators appearing
    # in any U_i), likewise V and W.  Then sum_i U'_i V'_i W'_i = L_U L_V L_W * (Brent tensor),
    # so the target is scaled by the same factor.  Python ints are far faster than Fractions.
    scale = 1
    if not algebraic_sym:
        from math import lcm
        Ls = []
        for k in range(3):
            L = 1
            for t in triads:
                for row in t[k]:
                    for v in row:
                        if v.denominator != 1: L = lcm(L, v.denominator)
            Ls.append(L)
        maxabs = max((abs(v) for t in triads for M in t for row in M for v in row), default=0)
        triads = [[[[int(v * Ls[k]) for v in row] for row in t[k]] for k in range(3)] for t in triads]
        scale = Ls[0] * Ls[1] * Ls[2]
        res['denominator_lcm'] = [str(L) for L in Ls]
    else:
        maxabs = 0
    # Fast path (exact): integer coefficients that provably fit in int64 throughout.
    # Each term index (a,b,c,d,e,f) is packed into one integer, products are formed per triad
    # by broadcasting, and equal indices are summed after a sort.  The bound below guarantees
    # that no product and no partial sum can exceed 2^62 in absolute value.
    fast = False
    if not algebraic_sym:
        import numpy as np
        mU = max((abs(v) for t in triads for row in t[0] for v in row), default=0)
        mV = max((abs(v) for t in triads for row in t[1] for v in row), default=0)
        mW = max((abs(v) for t in triads for row in t[2] for v in row), default=0)
        terms = sum(len(nz(U)) * len(nz(V)) * len(nz(W)) for U, V, W in triads)
        # Rigorous bound on every product and partial sum: sum over triads of (sum|u|)(sum|v|)(sum|w|)
        # (each product u*v*w is one term of that expansion, and any partial sum is at most the total).
        l1 = lambda M: sum(abs(v) for row in M for v in row)
        bound = sum(l1(U) * l1(V) * l1(W) for U, V, W in triads)
        res['terms'] = terms
        if bound + scale < 2 ** 62 and not os.environ.get('FMM_FORCE_CRT'):
            fast = True; moduli = [None]
        else:
            # Exact by the Chinese remainder theorem: every coefficient S of the accumulated tensor
            # and its target T satisfy |S - T| <= bound + scale.  The computation is repeated modulo
            # primes p < 2^30 (all products and partial sums then stay below 2^62) until the product
            # of the primes exceeds bound + scale.  S == T modulo every prime then forces S == T,
            # and S != T modulo any one prime proves S != T.
            fast = True; moduli = []; M = 1; q = 2 ** 30
            def isprime(n):
                if n % 2 == 0: return False
                i = 3
                while i * i <= n:
                    if n % i == 0: return False
                    i += 2
                return True
            while M <= bound + scale:
                q -= 1
                while not isprime(q): q -= 1
                moduli.append(q); M *= q
        res['int64_bound_log2'] = max(bound, 1).bit_length()
    dense = fast and (terms > 2 * 10 ** 8 or bool(os.environ.get('FMM_FORCE_DENSE')))
    if dense:
        # Dense exact accumulation for schemes with billions of terms.  The coefficient tensor is
        # indexed by (ab, cd, ef) with ab = a*Y+b, cd = c*Z+d, ef = e*Z+f, exactly as the packed key of
        # the sparse path.  It is processed in blocks of ab so that a block holds at most about 6e7
        # entries.  Each triad adds the outer product of its nonzero values into the block positions
        # np.ix_(iu, iv, iw); the index lists within one triad are distinct, so the in-place add is
        # exact.  A boolean array records which coefficients received at least one term, so wrong and
        # missing entries are counted with the same meaning as in the sparse path.
        YZ, ZX = Y * Z, Z * X
        blk = max(0, int(os.environ.get('FMM_DENSE_BLOCK', 0) or 0)) or max(1, min(X * Y, (6 * 10 ** 7) // (YZ * ZX)))
        sp = []
        for U, V, W in triads:
            nu, nv, nw = nz(U), nz(V), nz(W)
            if not (nu and nv and nw): continue
            sp.append((np.array([a * Y + b for a, b, _ in nu], dtype=np.int64), [u for _, _, u in nu],
                       np.array([c * Z + d for c, d, _ in nv], dtype=np.int64), [v for _, _, v in nv],
                       np.array([e * Z + f for f, e, _ in nw], dtype=np.int64), [w for _, _, w in nw]))
        bad = 0; missing = 0; touched_total = 0; worst = None
        for lo in range(0, X * Y, blk):
            hi = min(X * Y, lo + blk)
            touched = np.zeros((hi - lo, YZ, ZX), dtype=bool)
            wrong = np.zeros((hi - lo, YZ, ZX), dtype=bool)
            # target positions inside this block: ab=(a,b), cd=(b,d), ef=(e=a, f=d)
            tab = np.arange(lo, hi); ta, tb = tab // Y, tab % Y
            for k, mod in enumerate(moduli):
                T = np.zeros((hi - lo, YZ, ZX), dtype=np.int64)
                for ia, ua, ib, vb, ic, wc in sp:
                    sel = (ia >= lo) & (ia < hi)
                    if not sel.any(): continue
                    iu = ia[sel] - lo
                    if mod:
                        u = np.array([x % mod for x, s_ in zip(ua, sel) if s_], dtype=np.int64)
                        v = np.array([x % mod for x in vb], dtype=np.int64); w = np.array([x % mod for x in wc], dtype=np.int64)
                        o = ((u[:, None, None] * v[None, :, None]) % mod) * w[None, None, :] % mod
                        ix = np.ix_(iu, ib, ic)
                        T[ix] = (T[ix] + o) % mod
                    else:
                        u = np.array([x for x, s_ in zip(ua, sel) if s_], dtype=np.int64)
                        v = np.array(vb, dtype=np.int64); w = np.array(wc, dtype=np.int64)
                        T[np.ix_(iu, ib, ic)] += u[:, None, None] * v[None, :, None] * w[None, None, :]
                    if k == 0: touched[np.ix_(iu, ib, ic)] = True
                tgt_val = (scale % mod) if mod else scale
                tgt = np.zeros_like(T)
                for d_ in range(Z):
                    tgt[tab - lo, tb * Z + d_, ta * Z + d_] = tgt_val
                w_m = (T != tgt) & touched
                if worst is None and w_m.any():
                    q = np.argwhere(w_m)[0]; ab_, cd_, ef_ = int(q[0]) + lo, int(q[1]), int(q[2])
                    worst = ((ab_ // Y, ab_ % Y, cd_ // Z, cd_ % Z, ef_ // Z, ef_ % Z), ('%d mod %d' % (int(T[tuple(q)]), mod)) if mod else str(int(T[tuple(q)])), int(tgt[tuple(q)]))
                wrong |= w_m
                del T
            istgt = np.zeros((hi - lo, YZ, ZX), dtype=bool)
            for d_ in range(Z): istgt[tab - lo, tb * Z + d_, ta * Z + d_] = True
            missing += int((istgt & ~touched).sum()); bad += int(wrong.sum()); touched_total += int(touched.sum())
            del touched, wrong, istgt
        T = range(touched_total)
        res['path'] = 'dense_int64' if moduli == [None] else 'dense_crt'
        res['dense_block'] = blk
        if moduli != [None]: res['crt_moduli'] = moduli
    elif fast:
        CHUNK = 2 * 10 ** 7
        def reduce(ks, vs, mod):
            K = np.concatenate(ks); Vv = np.concatenate(vs)
            order = np.argsort(K, kind='stable'); K = K[order]; Vv = Vv[order]
            uk, start = np.unique(K, return_index=True)
            S = np.add.reduceat(Vv, start)
            return uk, (S % mod if mod else S)
        wrong = None; uk = None; worst = None
        for mod in moduli:
            red = (lambda x: x % mod) if mod else (lambda x: x)
            acc_k = np.zeros(0, dtype=np.int64); acc_v = np.zeros(0, dtype=np.int64)
            keys = []; vals = []; pending = 0
            for U, V, W in triads:
                nu, nv, nw = nz(U), nz(V), nz(W)
                if not (nu and nv and nw): continue
                ka = np.array([(a * Y + b) for a, b, _ in nu], dtype=np.int64); va = np.array([red(u) for _, _, u in nu], dtype=np.int64)
                kb = np.array([(c * Z + d) for c, d, _ in nv], dtype=np.int64); vb = np.array([red(v) for _, _, v in nv], dtype=np.int64)
                kc = np.array([(e * Z + f) for f, e, _ in nw], dtype=np.int64); vc = np.array([red(w) for _, _, w in nw], dtype=np.int64)
                # index = ((a*Y+b) * (Y*Z) + (c*Z+d)) * (X*Z) + (e*Z+f)
                k = ((ka[:, None, None] * (Y * Z) + kb[None, :, None]) * (X * Z) + kc[None, None, :]).ravel()
                if mod:
                    v = (((va[:, None, None] * vb[None, :, None]) % mod) * vc[None, None, :] % mod).ravel()
                else:
                    v = (va[:, None, None] * vb[None, :, None] * vc[None, None, :]).ravel()
                keys.append(k); vals.append(v); pending += len(k)
                if pending >= CHUNK:
                    acc_k, acc_v = reduce([acc_k] + keys, [acc_v] + vals, mod); keys = []; vals = []; pending = 0
            uk_m, sums = reduce([acc_k] + keys, [acc_v] + vals, mod) if (keys or len(acc_k)) else (np.zeros(0, dtype=np.int64), np.zeros(0, dtype=np.int64))
            if uk is None:
                uk = uk_m
                ef = uk % (X * Z); rest = uk // (X * Z); cd = rest % (Y * Z); ab = rest // (Y * Z)
                a = ab // Y; b = ab % Y; c = cd // Z; d = cd % Z; e = ef // Z; f = ef % Z
                is_target = (b == c) & (a == e) & (d == f)
            elif not np.array_equal(uk, uk_m):
                raise RuntimeError('key sets differ between moduli')
            tgt = np.where(is_target, (scale % mod) if mod else scale, 0)
            w_m = sums != tgt
            if worst is None and w_m.any():
                w0 = int(np.nonzero(w_m)[0][0])
                worst = ((int(a[w0]), int(b[w0]), int(c[w0]), int(d[w0]), int(e[w0]), int(f[w0])), ('%d mod %d' % (int(sums[w0]), mod)) if mod else str(int(sums[w0])), int(tgt[w0]))
            wrong = w_m if wrong is None else (wrong | w_m)
        bad = int(wrong.sum()) if wrong is not None else 0
        missing = X * Y * Z - int(is_target.sum()) if uk is not None else X * Y * Z
        T = uk if uk is not None else []
        res['path'] = 'numpy_int64' if moduli == [None] else 'numpy_crt'
        if moduli != [None]: res['crt_moduli'] = moduli
    else:
        T = defaultdict(int)
        for U, V, W in triads:
            nu, nv, nw = nz(U), nz(V), nz(W)
            for a, b, u in nu:
                for c, d, v in nv:
                    uv = u * v
                    for f, e, w in nw:
                        T[(a, b, c, d, e, f)] += uv * w
        # compare with the matrix multiplication tensor
        bad = 0; worst = None; checked = 0
        for key, val in T.items():
            a, b, c, d, e, f = key
            target = scale if (b == c and a == e and d == f) else 0
            if algebraic_sym:
                import sympy
                val = sympy.simplify(val - target); target = 0
            if val != target:
                bad += 1
                if worst is None: worst = (key, str(val), target)
        # entries of the target that T never touched (must be 1) -> missing
        missing = 0
        for a in range(X):
            for b in range(Y):
                for d in range(Z):
                    if (a, b, b, d, a, d) not in T: missing += 1
        res['path'] = 'python_exact'
    res['nonzero_terms'] = len(T); res['max_abs_coeff'] = str(maxabs)
    if bad or missing:
        res['verdict'] = 'FAIL'; res['wrong_entries'] = bad; res['missing_entries'] = missing; res['first_wrong'] = worst
    else:
        res['verdict'] = 'PASS' if len(triads) == stated_rank else 'PASS_RANK_MISMATCH'
    return res

if __name__ == "__main__":
    path = sys.argv[1]; X, Y, Z, r = map(int, sys.argv[2:6])
    print(json.dumps(verify(path, X, Y, Z, r)))
