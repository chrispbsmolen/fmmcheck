"""Lean exact check of a https://github.com/khoruzhii/lita scheme read directly from its .npz (CSR rows of U, V, W
over row-major flattened N x N matrices, C = A.B with C_k = sum_q U[q,i] V[q,j] W[q,k] A_i B_j).  Used because
converting the largest of these dense schemes to text needs more memory than verify_fmm.py has here.

Method: every factor is scaled by the lcm of its denominators, so all coefficients are integers.  The tensor
sum_q U[q,i] V[q,j] W[q,k] is accumulated densely, in blocks of i, modulo primes p < 2^30 until the product of
the primes exceeds bound + scale (bound = sum_q |U_q|_1 |V_q|_1 |W_q|_1), and compared with scale * T_ijk,
T_ijk = [i=(a,b), j=(b,d), k=(a,d)].  Equal modulo every prime forces equality; a difference modulo one prime
proves a wrong coefficient.  Cross-checked against verify_fmm.py in tools/lita_exact_test.py.
    python3 tools/lita_exact.py SCHEME.npz [--corrupt]"""
import sys, json, math
import numpy as np
path = sys.argv[1]; z = np.load(path, allow_pickle=False); meta = json.loads(z['metadata_json'].item())
N = meta['tensor'][0]; R = meta['rank']; N2 = N * N
def axis(name):
    ip, ix = z[name + '_indptr'], z[name + '_indices'].astype(np.int64)
    num, den = [int(x) for x in z[name + '_numerators']], [int(x) for x in z[name + '_denominators']]
    L = 1
    for d in set(den): L = L * d // math.gcd(L, d)
    vals = [n * (L // d) for n, d in zip(num, den)]
    return ip, ix, vals, L
(uip, uix, uval, Lu), (vip, vix, vval, Lv), (wip, wix, wval, Lw) = axis('u'), axis('v'), axis('w')
if '--corrupt' in sys.argv: uval[0] += Lu          # change one coefficient of U_0 by 1 (after scaling)
scale = Lu * Lv * Lw
l1 = lambda ip, vals, q: sum(abs(x) for x in vals[ip[q]:ip[q + 1]])
bound = sum(l1(uip, uval, q) * l1(vip, vval, q) * l1(wip, wval, q) for q in range(R))
moduli = []; M = 1; p = 2 ** 30
def isprime(n):
    if n % 2 == 0: return False
    i = 3
    while i * i <= n:
        if n % i == 0: return False
        i += 2
    return True
while M <= bound + scale:
    p -= 1
    while not isprime(p): p -= 1
    moduli.append(p); M *= p
import os
blk = max(0, int(os.environ.get('LITA_BLOCK', 0) or 0)) or max(1, min(N2, (4 * 10 ** 7) // (N2 * N2)))
prog = path + ('.corrupt' if '--corrupt' in sys.argv else '') + '.progress.json'   # resumable: one entry per finished block
import os
done = json.load(open(prog)) if os.path.exists(prog) and json.load(open(prog)).get('blk') == blk else {'blk': blk, 'blocks': {}}
wrong = 0; missing = 0; touched_total = 0
for lo in range(0, N2, blk):
    if str(lo) in done['blocks']:
        b = done['blocks'][str(lo)]; wrong += b[0]; missing += b[1]; touched_total += b[2]; continue
    hi = min(N2, lo + blk); touched = np.zeros((hi - lo, N2, N2), dtype=bool); bad = np.zeros_like(touched)
    ia = np.arange(lo, hi); ta, tb = ia // N, ia % N
    for k, mod in enumerate(moduli):
        T = np.zeros((hi - lo, N2, N2), dtype=np.int64)
        for q in range(R):
            iu = uix[uip[q]:uip[q + 1]]; sel = (iu >= lo) & (iu < hi)
            if not sel.any(): continue
            u = np.array([x % mod for x, s in zip(uval[uip[q]:uip[q + 1]], sel) if s], dtype=np.int64)
            v = np.array([x % mod for x in vval[vip[q]:vip[q + 1]]], dtype=np.int64)
            w = np.array([x % mod for x in wval[wip[q]:wip[q + 1]]], dtype=np.int64)
            ix = np.ix_(iu[sel] - lo, vix[vip[q]:vip[q + 1]], wix[wip[q]:wip[q + 1]])
            T[ix] = (T[ix] + ((u[:, None, None] * v[None, :, None]) % mod) * w[None, None, :]) % mod
            if k == 0: touched[ix] = True
        tv = scale % mod
        for d in range(N):                        # target positions i=(a,b), j=(b,d), k=(a,d)
            T[ia - lo, tb * N + d, ta * N + d] -= tv
        T %= mod
        bad |= (T != 0) & touched
        del T
    ist = np.zeros_like(touched)
    for d in range(N): ist[ia - lo, tb * N + d, ta * N + d] = True
    bm, bw, bt = int((ist & ~touched).sum()), int(bad.sum()), int(touched.sum())
    missing += bm; wrong += bw; touched_total += bt
    done['blocks'][str(lo)] = [bw, bm, bt]; json.dump(done, open(prog, 'w'))
    print('block', lo, 'of', N2, bw, bm, file=sys.stderr, flush=True)
    del touched, bad, ist
verdict = 'FAIL' if (wrong or missing) else 'PASS'
print(json.dumps(dict(scheme=path.split('/')[-1], N=N, rank=R, verdict=verdict, wrong_entries=wrong, missing_entries=missing,
                      nonzero_coefficients=touched_total, crt_moduli=len(moduli), corrupted='--corrupt' in sys.argv)))
