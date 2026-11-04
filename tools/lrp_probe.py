"""Floating point screen of a catalogue _LRP.mpl file (L: r x XY, R: r x YZ, P: XZ x r) over the 8 row/column-major
conventions for vec(A), vec(B), vec(C): prints the largest error against A.B for random integer A, B.  Used only to
find the convention; the converted scheme is then verified exactly.  Run in the folder holding FORMAT_LRP.mpl.bz2:
    python3 lrp_probe.py FORMAT"""
import re, bz2, random, itertools
from fractions import Fraction
import numpy as np
import sys
f = sys.argv[1]; X, Y, Z = map(int, f.split('x'))
s = bz2.open(f + '_LRP.mpl.bz2', 'rt').read()
mats = []
for m in re.finditer(r'Matrix\((\d+), (\d+), \[\[(.*?)\]\]\)', s, re.S):
    rows = m.group(3).replace('\n', '').split('],[')
    mats.append(np.array([[float(Fraction(x)) for x in row.split(',')] for row in rows]))
L, R, P = mats; print(f, 'shapes', L.shape, R.shape, P.shape)
rng = np.random.default_rng(1); A = rng.integers(-5, 6, (X, Y)).astype(float); B = rng.integers(-5, 6, (Y, Z)).astype(float); C = A @ B
def vec(M, o): return M.reshape(-1) if o == 'r' else M.T.reshape(-1)
for oa, ob in itertools.product('rc', repeat=2):
    m = (L @ vec(A, oa)) * (R @ vec(B, ob)); c = P @ m
    for oc in 'rc':
        for tgt, name in ((C, 'AB'),):
            err = np.abs(c - vec(tgt, oc)).max()
            print(f, oa, ob, oc, name, '%.3g' % err)
    # also try L acting on A^T etc via other shapes
