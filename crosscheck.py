"""Separately written second evaluator: its own parser (character scanner, nothing shared with verify_fmm.py) and
a different method.  It evaluates the file's identity A.B = sum_i Tr(U_i^T A) Tr(V_i^T B) W_i^T for random A, B
modulo the prime p = 2^61 - 1 (rational entries are reduced mod p; free parameters get random values mod p) and
counts the wrong entries of A.B.  A wrong entry mod p proves the file wrong (if the exact identity held, it would
hold mod p).  All entries right on random inputs is strong evidence, not a proof; verify_fmm.py is the proof.
    python3 crosscheck.py FILE X Y Z [trials]          memory: about 1 GB for the largest catalogue files"""
import sys, bz2, random, re, ast
p = (1 << 61) - 1
path = sys.argv[1]; X, Y, Z = map(int, sys.argv[2:5]); trials = int(sys.argv[5]) if len(sys.argv) > 5 else 2
txt = bz2.open(path, 'rt').read() if path.endswith('.bz2') else open(path).read()
txt = txt.replace('\\\n', '')
start = txt.index('TriadSet(')
rng = random.Random(61); params = {}
def value(tok):
    tok = ''.join(tok.split())
    if re.fullmatch(r'-?\d+', tok): return int(tok) % p
    m = re.fullmatch(r'(-?\d+)/(\d+)', tok)
    if m: return int(m.group(1)) * pow(int(m.group(2)), p - 2, p) % p
    def ev(n):   # + - * / ^ on integers and names, mod p
        if isinstance(n, ast.Expression): return ev(n.body)
        if isinstance(n, ast.Constant): return int(n.value) % p
        if isinstance(n, ast.Name):
            if n.id not in params: params[n.id] = rng.randrange(2, p - 1)
            return params[n.id]
        if isinstance(n, ast.UnaryOp): return (-ev(n.operand)) % p if isinstance(n.op, ast.USub) else ev(n.operand)
        if isinstance(n, ast.BinOp):
            a, b = ev(n.left), ev(n.right)
            if isinstance(n.op, ast.Add): return (a + b) % p
            if isinstance(n.op, ast.Sub): return (a - b) % p
            if isinstance(n.op, ast.Mult): return a * b % p
            if isinstance(n.op, ast.Div): return a * pow(b, p - 2, p) % p
            if isinstance(n.op, (ast.Pow, ast.BitXor)): return pow(a, int(ast.literal_eval(n.right)), p)
        raise ValueError('unsupported ' + type(n).__name__)
    return ev(ast.parse(tok.replace('^', '**'), mode='eval'))
def matrices():
    i = start
    while True:
        j = txt.find('Matrix(', i)
        if j < 0: return
        k = txt.index('[', j); r, c = (int(x) for x in txt[j + 7:k].split(',')[:2])
        depth = 0; q = k
        while True:                                   # scan to the bracket that closes the row list
            ch = txt[q]
            if ch == '[': depth += 1
            elif ch == ']':
                depth -= 1
                if depth == 0: break
            q += 1
        rows = [row for row in re.split(r'\]\s*,\s*\[', txt[k + 1:q].strip()[1:-1])]
        M = [[value(t) for t in row.split(',')] for row in rows]
        if len(M) != r or any(len(x) != c for x in M): raise ValueError('shape at %d' % j)
        yield M; i = q
mats = list(matrices()); assert len(mats) % 3 == 0, 'matrix count %d' % len(mats)
triads = [mats[3 * t: 3 * t + 3] for t in range(len(mats) // 3)]
for t in range(trials):
    A = [[rng.randrange(p) for _ in range(Y)] for _ in range(X)]; B = [[rng.randrange(p) for _ in range(Z)] for _ in range(Y)]
    C = [[0] * Z for _ in range(X)]
    for U, V, W in triads:
        assert len(U) == X and len(U[0]) == Y and len(V) == Y and len(V[0]) == Z and len(W) == Z and len(W[0]) == X
        s = sum(U[a][b] * A[a][b] for a in range(X) for b in range(Y)) % p
        if not s: continue
        r = sum(V[b][d] * B[b][d] for b in range(Y) for d in range(Z)) % p
        if not r: continue
        sr = s * r % p
        for d in range(Z):
            Wd = W[d]
            for a in range(X):
                if Wd[a]: C[a][d] = (C[a][d] + sr * Wd[a]) % p
    wrong = sum(1 for a in range(X) for d in range(Z) if C[a][d] != sum(A[a][b] * B[b][d] for b in range(Y)) % p)
    print('%s triads %d trial %d: %d of %d entries of A.B wrong mod 2^61-1' % (path.split('/')[-1], len(triads), t, wrong, X * Z), flush=True)
