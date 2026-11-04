"""Write results/ledger.csv (one row per format), results/rank_mismatches.csv, results/defects.csv and
results/summary.json from data/entries.json, data/page_ranks.json, data/results.jsonl (this run, current
verify_fmm.py) and, for formats this run did not reach, data/results_previous.jsonl results/results_previous.jsonl (earlier runs, 2-3 Nov 2026,
with earlier versions of verify_fmm.py that are not in this repository, all exact; the 'run' column says which).
data/results.jsonl is copied to results/results_this_run.jsonl.  Where both exist they are compared.
Every number in the README and in results/VERIFICATION.md is read from summary.json."""
import json, csv, collections, os, bz2, hashlib, subprocess
W = os.path.dirname(os.path.abspath(__file__)); D = W + '/data'; O = W + '/results'; os.makedirs(O, exist_ok=True)
E = json.load(open(D + '/entries.json'))
P = json.load(open(D + '/page_ranks.json')) if os.path.exists(D + '/page_ranks.json') else {}
R = {}; PREV = {}
THIS = D + '/results.jsonl' if os.path.exists(D + '/results.jsonl') else O + '/results_this_run.jsonl'
for l in open(THIS):
    r = json.loads(l); r['run'] = 'this run'; R[r['fmt']] = r
PREVF = O + '/results_previous.jsonl'
if os.path.exists(PREVF):
    for l in open(PREVF):
        r = json.loads(l); PREV[r['fmt']] = r
key = lambda r: (r['verdict'], r.get('rank'), r.get('wrong_entries'), r.get('missing_entries'), r.get('file_sha256'))
overlap = [f for f in R if f in PREV]; disagree = [f for f in overlap if key(R[f]) != key(PREV[f])]
for f, r in PREV.items():
    if f not in R: R[f] = r
def no_scheme_kind(f):
    t = bz2.open(D + '/tensors/%s_tensor.mpl.bz2' % f, 'rt').read()
    if 'Tensor:=Tensor:' in t: return 'placeholder (Tensor:=Tensor:, no triads)'
    if 'ReturnFastMatrixMultiplication' in t and 'TriadSet' not in t: return 'stub (calls ReturnFastMatrixMultiplication, no triads)'
    if 'TriadSet' in t and not t.rstrip().endswith(';'): return 'truncated file (text ends inside a matrix)'
    return 'unreadable: ' + R[f].get('note', '')
rows = []
for f in sorted(E, key=lambda f: (E[f]['X'], E[f]['Y'], E[f]['Z'])):
    r = R.get(f, {}); pr = P.get(f, {}).get('page_rank'); fr = r.get('rank'); v = r.get('verdict', 'NOT_RUN')
    cmp_ = ''
    if v == 'PASS_RANK_MISMATCH':
        cmp_ = 'page = file' if pr == fr else ('page = front' if pr == E[f]['rank'] else ('format page empty or gives no rank' if pr is None else 'page differs from both'))
    detail = no_scheme_kind(f) if v == 'UNPARSED' else ('%d wrong, %d missing' % (r.get('wrong_entries', 0), r.get('missing_entries', 0)) if v == 'FAIL' else '')
    rows.append(dict(format=f, front_page_rank=E[f]['rank'], format_page_rank='' if pr is None else pr, file_triads='' if fr is None else fr,
                     verdict=v, detail=detail, rank_comparison=cmp_, method=r.get('path') or ('none (no scheme)' if v == 'UNPARSED' else ''), verifier=r.get('verifier_sha256') or 'earlier version',
                     run=r.get('run', ''), file_sha256=r.get('file_sha256', ''), front_page_source=E[f].get('source', '')))
fields = list(rows[0])
for name, sel in (('ledger.csv', rows), ('rank_mismatches.csv', [x for x in rows if x['verdict'] == 'PASS_RANK_MISMATCH']),
                  ('defects.csv', [x for x in rows if x['verdict'] in ('FAIL', 'UNPARSED', 'BADSHAPE', 'ERROR', 'NOT_RUN')])):
    with open(O + '/' + name, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields); w.writeheader(); w.writerows(sel)
mm = [x for x in rows if x['verdict'] == 'PASS_RANK_MISMATCH']
try: commit = subprocess.run(['git', 'log', '-1', '--format=%h', '--', 'verify_fmm.py'], cwd=W, capture_output=True, text=True).stdout.strip()
except Exception: commit = ''
S = dict(snapshot=open(D + '/SNAPSHOT_DATE').read().strip() if os.path.exists(D + '/SNAPSHOT_DATE') else '',
         formats=len(E), run=len(R), parametric_this_run=sum(1 for f in R if R[f].get('run') == 'this run' and R[f].get('ring') == 'parametric'), runs=dict(collections.Counter(x['run'] for x in rows)), rerun_overlap=len(overlap), rerun_disagreements=disagree, not_run=sum(1 for x in rows if x['verdict'] == 'NOT_RUN'),
         verdicts=dict(collections.Counter(x['verdict'] for x in rows)),
         mismatch=dict(collections.Counter(x['rank_comparison'] for x in mm)),
         file_fewer_triads_than_front=sum(1 for x in mm if x['file_triads'] < x['front_page_rank']),
         file_more_triads_than_front=sum(1 for x in mm if x['file_triads'] > x['front_page_rank']),
         no_scheme=dict(collections.Counter(x['detail'] for x in rows if x['verdict'] == 'UNPARSED')),
         fails=[{k: x[k] for k in ('format', 'front_page_rank', 'format_page_rank', 'file_triads', 'detail', 'front_page_source', 'file_sha256')} for x in rows if x['verdict'] == 'FAIL'],
         methods=dict(collections.Counter(x['method'] for x in rows)),
         verifiers_this_run=dict(collections.Counter(R[f].get('verifier_sha256') for f in R if R[f].get('run') == 'this run')), verifier_sha256_now=hashlib.sha256(open(W + '/verify_fmm.py', 'rb').read()).hexdigest()[:16], verifier_commit=commit,
         seconds_total=round(sum(R[f].get('seconds', 0) for f in R if R[f].get('run') == 'this run')))
json.dump(S, open(O + '/summary.json', 'w'), indent=1)
if THIS != O + '/results_this_run.jsonl':
    import shutil; shutil.copyfile(THIS, O + '/results_this_run.jsonl')
print(json.dumps({k: S[k] for k in ('formats', 'run', 'verdicts', 'mismatch', 'no_scheme')}, indent=1))
