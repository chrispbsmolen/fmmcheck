"""Cross-check tools/lita_exact.py against verify_fmm.py (via tools/lita_convert.py) on the smallest lita scheme,
original and with one coefficient changed."""
import sys, os, json, subprocess
W = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'); sys.path.insert(0, W)
from verify_fmm import verify
npz = sys.argv[1]; N = int(os.path.basename(npz).split('x')[0]); r = int(npz.split('_r')[1].split('.')[0])
a = json.loads(subprocess.run([sys.executable, W + '/tools/lita_exact.py', npz], capture_output=True, text=True).stdout)
b = json.loads(subprocess.run([sys.executable, W + '/tools/lita_exact.py', npz, '--corrupt'], capture_output=True, text=True).stdout)
tmp = W + '/data/lita_test_%d.mpl.bz2' % N
subprocess.run([sys.executable, W + '/tools/lita_convert.py', npz, tmp], check=True, capture_output=True)
v = verify(tmp, N, N, N, r)
subprocess.run([sys.executable, W + '/tools/lita_convert.py', npz, tmp + '.c.mpl.bz2', '--corrupt'], check=True, capture_output=True)
vc = verify(tmp + '.c.mpl.bz2', N, N, N, r)
print('scheme:', os.path.basename(npz))
print('lita_exact original :', a['verdict'], a['wrong_entries'], a['missing_entries'], a['nonzero_coefficients'])
print('verify_fmm original :', v['verdict'], v.get('wrong_entries', 0), v.get('missing_entries', 0), v.get('nonzero_terms'))
print('lita_exact corrupted:', b['verdict'], b['wrong_entries'], b['missing_entries'])
print('verify_fmm corrupted:', vc['verdict'], vc.get('wrong_entries', 0), vc.get('missing_entries', 0))
ok = a['verdict'] == ('PASS' if v['verdict'] in ('PASS', 'PASS_RANK_MISMATCH') else 'FAIL') and a['nonzero_coefficients'] == v.get('nonzero_terms') and b['verdict'] == 'FAIL' and (b['wrong_entries'], b['missing_entries']) == (vc.get('wrong_entries'), vc.get('missing_entries'))
print('agree' if ok else 'PROBLEM'); sys.exit(0 if ok else 1)
