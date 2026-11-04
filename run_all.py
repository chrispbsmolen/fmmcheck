"""Verify every downloaded scheme file.  python3 run_all.py [workers] [FORMATS_FILE]   (default 1 worker)
With FORMATS_FILE (one format per line) only those formats are checked.
Appends one JSON line per format to data/results.jsonl (resumable).  Formats with X*Y*Z below 8000 run on
the given number of workers; larger ones run one at a time.  Memory: up to about 4 GB per worker for formats
near X*Y*Z = 8000 and 5 GB for the largest, so use 1 or 2 workers on a machine with 8 GB."""
import json, os, sys, time, hashlib, multiprocessing as mp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
D = os.path.join(HERE, 'data'); E = json.load(open(D + '/entries.json')); OUT = D + '/results.jsonl'
VSHA = hashlib.sha256(open(os.path.join(HERE, 'verify_fmm.py'), 'rb').read()).hexdigest()[:16]
def work(fmt):
    from verify_fmm import verify
    e = E[fmt]; p = D + '/tensors/%s_tensor.mpl.bz2' % fmt; t0 = time.time()
    try: r = verify(p, e['X'], e['Y'], e['Z'], e['rank'])
    except Exception as ex: r = dict(verdict='ERROR', note=repr(ex)[:300])
    r['fmt'] = fmt; r['seconds'] = round(time.time() - t0, 2)
    r['file_sha256'] = hashlib.sha256(open(p, 'rb').read()).hexdigest()
    r['verifier_sha256'] = VSHA
    return r
if __name__ == '__main__':
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT): done.add(json.loads(l)['fmt'])
    todo = sorted((f for f in E if f not in done and os.path.exists(D + '/tensors/%s_tensor.mpl.bz2' % f)), key=lambda f: E[f]['X'] * E[f]['Y'] * E[f]['Z'])
    nw = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    if len(sys.argv) > 2:
        only = set(open(sys.argv[2]).read().split()); todo = [f for f in todo if f in only]
    print('to verify:', len(todo), flush=True)
    small = [f for f in todo if E[f]['X'] * E[f]['Y'] * E[f]['Z'] < 8000]; large = [f for f in todo if f not in small]
    k = 0
    with open(OUT, 'a') as out:
        for batch, n in ((small, nw), (large, 1)):
            with mp.Pool(n, maxtasksperchild=10) as pool:
                for r in pool.imap_unordered(work, batch):
                    out.write(json.dumps(r) + '\n'); out.flush(); k += 1
                    if r['verdict'] not in ('PASS', 'PASS_RANK_MISMATCH'): print('  ', r['fmt'], r['verdict'], r.get('note', ''), flush=True)
                    if k % 100 == 0: print(time.strftime('%H:%M'), k, '/', len(todo), flush=True)
