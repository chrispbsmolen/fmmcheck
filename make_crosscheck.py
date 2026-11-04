"""Run crosscheck.py on the six FAIL files, one PASS control of similar size, and the three files rebuilt by
make_fixes.py; write results/crosscheck_output.txt.  One at a time (the largest needs about 3 GB)."""
import subprocess, sys, os, json
W = os.path.dirname(os.path.abspath(__file__)); D = W + '/data'
S = json.load(open(W + '/results/summary.json'))
jobs = [(D + '/tensors/%s_tensor.mpl.bz2' % f['format'], f['format']) for f in S['fails']] + [(D + '/tensors/12x28x28_tensor.mpl.bz2', '12x28x28')] + \
       [(D + '/fixes/%s_tensor_rebuilt.mpl.bz2' % f, f) for f in ('11x25x28', '15x20x26', '8x27x30')]
out = []
for path, fmt in jobs:
    X, Y, Z = fmt.split('x')
    r = subprocess.run([sys.executable, W + '/crosscheck.py', path, X, Y, Z, '1'], capture_output=True, text=True)
    line = (r.stdout.strip() or ('ERROR ' + r.stderr.strip()[-200:])).replace(path.split('/')[-1], os.path.relpath(path, W))
    out.append(line); print(line, flush=True)
open(W + '/results/crosscheck_output.txt', 'w').write('\n'.join(out) + '\n')
