"""Evaluate a catalogue _raw.mpl file (MUL: m_k = (lin A)*(lin B); ADD: C_i_j = sum of m_k) on random
rational A, B and random parameter values, exactly, and compare with A.B.  Expressions are parsed with
Python's ast after replacing names; only + - * / and numbers and names are allowed."""
import sys, bz2, re, ast, random, operator
from fractions import Fraction
f = sys.argv[1]; X, Y, Z = map(int, f.split('x'))
s = bz2.open(f + '_raw.mpl.bz2', 'rt').read().replace('\n', '')
mul = s[s.find('MUL:=[') + 6: s.find(']:', s.find('MUL:=['))]
add = s[s.find('ADD:=[') + 6: s.find(']:', s.find('ADD:=['))]
rng = random.Random(7)
env = {}
for i in range(X):
    for j in range(Y): env['A_%d_%d' % (i + 1, j + 1)] = Fraction(rng.randint(-9, 9))
for i in range(Y):
    for j in range(Z): env['B_%d_%d' % (i + 1, j + 1)] = Fraction(rng.randint(-9, 9))
for p in set(re.findall(r'\ba\d+\b', s)): env[p] = Fraction(rng.randint(2, 50), rng.randint(2, 50))
OPS = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul, ast.Div: operator.truediv}
def ev(node):
    if isinstance(node, ast.Expression): return ev(node.body)
    if isinstance(node, ast.BinOp): return OPS[type(node.op)](ev(node.left), ev(node.right))
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub): return -ev(node.operand)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.UAdd): return ev(node.operand)
    if isinstance(node, ast.Constant): return Fraction(node.value)
    if isinstance(node, ast.Name): return env[node.id]
    raise ValueError(type(node).__name__)
def split_eqs(txt):
    out = []; depth = 0; cur = ''
    for ch in txt:
        if ch == '(': depth += 1
        if ch == ')': depth -= 1
        if ch == ',' and depth == 0: out.append(cur); cur = ''
        else: cur += ch
    if cur: out.append(cur)
    return out
for eq in split_eqs(mul):
    k, e = eq.split('=', 1); env[k.strip()] = ev(ast.parse(e, mode='eval'))
wrong = 0; n = 0
for eq in split_eqs(add):
    k, e = eq.split('=', 1); k = k.strip()
    if not k.startswith('C_'): continue
    i, j = map(int, k[2:].split('_')); n += 1
    true = sum(env['A_%d_%d' % (i, t)] * env['B_%d_%d' % (t, j)] for t in range(1, Y + 1))
    if ev(ast.parse(e, mode='eval')) != true: wrong += 1
print(f, 'products', len(split_eqs(mul)), 'C entries', n, 'wrong', wrong)
