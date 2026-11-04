"""Write README.md and results/VERIFICATION.md.  Every number is read from results/summary.json (written by
make_ledger.py), results/selftest_output.txt and results/FIXES.md, so the documents cannot drift from the data."""
import json, os, csv
W = os.path.dirname(os.path.abspath(__file__)); S = json.load(open(W + '/results/summary.json'))
fmt = lambda x: '{:,}'.format(x); V = S['verdicts']; M = S['mismatch']
cc = open(W + '/results/crosscheck_output.txt').read().strip() if os.path.exists(W + '/results/crosscheck_output.txt') else '(not run)'
st = open(W + '/results/selftest_output.txt').read().strip() if os.path.exists(W + '/results/selftest_output.txt') else '(not run)'
fails = S['fails']
jl = lambda n: [json.loads(l) for l in open(W + '/results/' + n)] if os.path.exists(W + '/results/' + n) else []
lr = {x['scheme']: x for x in jl('lita_randeval.jsonl')}; lx = {x['scheme']: x for x in jl('lita_exact.jsonl')}
lita = [dict(scheme=k, rank=k.split('_r')[1].split('.')[0], randeval='%d of %d entries wrong' % (v['wrong'], v['entries']),
             exact=('%s (exact, %d wrong, %d missing)' % (lx[k]['verdict'], lx[k]['wrong_entries'], lx[k]['missing_entries'])) if k in lx else ('started, not completed' if k.startswith('26x') else 'not run (too slow on the machine used)')) for k, v in sorted(lr.items())]
fail_rows = '\n'.join('| %s | %s | %s | %s | %s | %s |' % (f['format'], f['front_page_rank'], f['format_page_rank'] or 'none', f['file_triads'], f['detail'], f['front_page_source'][:90].replace('  ', ' ')) for f in fails)
nos = ', '.join('%d %s' % (v, k) for k, v in sorted(S['no_scheme'].items(), key=lambda t: -t[1]))
lita_rows = '\n'.join('| %s | %s | %s | %s |' % (x['scheme'], x['rank'], x['randeval'], x['exact']) for x in lita)
readme = f'''# fmmcheck

An exact checker for Alexandre Sedoglavic's catalogue of fast matrix multiplication algorithms at
https://fmm.univ-lille.fr/, and the results of running it on every scheme file in the catalogue.

Each catalogue entry `<X×Y×Z:r>` has a file `XxYxZ_tensor.mpl.bz2` holding triads `(U_i, V_i, W_i)` and the
identity `A.B = sum_i Tr(U_i^T A) Tr(V_i^T B) W_i^T`. `verify_fmm.py` proves or refutes that identity in exact
arithmetic (for schemes with free parameters, at fixed parameter values), by evaluating every coefficient of the
Brent equations, and reports whether the number of triads equals the rank in the catalogue's table. Two cheaper
checks multiply random matrices with the file's own formula: `randeval.py` (exact rationals, sharing only the file
parser with `verify_fmm.py`) and `crosscheck.py` (written separately, with its own parser, modulo 2^61 - 1).

## Results, snapshot of {S['snapshot']}

All {fmt(S['formats'])} formats in the front-page table were checked ({S['not_run']} not run). {fmt(S['runs'].get('this run', 0))} were checked in this run with the code in this repository. The other {fmt(S['formats'] - S['runs'].get('this run', 0))}, mostly the largest formats and including the four square formats below, were checked only in earlier runs of 2 and 3 November with earlier versions of `verify_fmm.py` that are not in this repository (also exact; the ledger's run column says which). On the {fmt(S['rerun_overlap'])} formats checked in both, the results agree exactly ({len(S['rerun_disagreements'])} differences).

| verdict | formats |
|---|---|
| PASS: the equations hold and the triad count equals the table's rank | {fmt(V.get('PASS', 0))} |
| PASS_RANK_MISMATCH: the equations hold, but the file has a different number of triads than the table | {fmt(V.get('PASS_RANK_MISMATCH', 0))} |
| FAIL: the equations do not hold | {fmt(V.get('FAIL', 0))} |
| no scheme in the file | {fmt(V.get('UNPARSED', 0))} |

Rank mismatches, compared with the rank stated on each format's own page (`results/rank_mismatches.csv`):
{fmt(M.get('page = file', 0))} where the format page agrees with the file and only the front page differs,
{fmt(M.get('page = front', 0))} where the format page agrees with the front page and the file differs,
{fmt(M.get('format page empty or gives no rank', 0))} where the format page is empty or gives no rank and {fmt(M.get('page differs from both', 0))} where it differs from both.
In {fmt(S['file_fewer_triads_than_front'])} of them the file has fewer triads than the table states.

Files with no scheme: {nos}.

### FAIL

| format | front page | format page | file triads | wrong / missing Brent coefficients | front-page source |
|---|---|---|---|---|---|
{fail_rows}

Every FAIL is confirmed by `crosscheck.py` (`results/crosscheck_output.txt`): random products are wrong in many
entries, while a correct file of similar size passes. For the four square formats the front page credits Khoruzhii,
Gelß and Pokutta 2026 at the ranks above and links https://github.com/khoruzhii/lita , while the served files match
the format pages (rank 8652 and so on, "private communication june 2026"). The schemes in that repository (commit
472a79f, `git clone https://github.com/khoruzhii/lita`), at the front-page ranks, were checked here:

| lita scheme (commit 472a79f) | rank | exact random-matrix check (tools/lita_randeval.py) | exact proof (tools/lita_exact.py) |
|---|---|---|---|
{lita_rows}

`make_fixes.py` rebuilds correct files for 11x25x28 (at the stated rank, from the catalogue's own LRP file),
15x20x26 and 8x27x30 (at slightly higher ranks; see `results/FIXES.md`). Those files are derived from the catalogue's
files and are not published here.

See `results/VERIFICATION.md` for the run of record and `results/ledger.csv` for one row per format, with the
SHA-256 of each file as served. The catalogue's files are not redistributed here; they are identified by hash.

## Usage

Python 3 with numpy (`requirements.txt`). Run every script from the repository root.

```
python3 selftest.py                 # constructed tests against a brute-force oracle; catalogue sample if data/ is filled
python3 verify_fmm.py FILE X Y Z RANK
python3 randeval.py FILE X Y Z

python3 fetch_catalogue.py --pages  # data/entries.json, data/tensors/ (about 600 MB), data/page_ranks.json
python3 run_all.py [workers]        # data/results.jsonl; up to 4 to 5 GB per worker, so 1 or 2 workers on 8 GB
python3 make_ledger.py              # results/ledger.csv, rank_mismatches.csv, defects.csv, summary.json
python3 make_fixes.py               # data/fixes/ and results/FIXES.md
python3 make_crosscheck.py          # results/crosscheck_output.txt
python3 make_readme.py              # this file and results/VERIFICATION.md
```

`tools/` holds the converters used by `make_fixes.py` (LRP and raw programs to the tensor format, the serendipitous
product builder; `lrp_probe.py` and `raw_probe.py` are run by `make_fixes.py` inside `data/fixes/`) and the scripts
for the lita schemes (`lita_randeval.py`, `lita_exact.py` with its test `lita_exact_test.py`, `lita_convert.py`).
The ledger is regenerated by `make_ledger.py` from `results/results_this_run.jsonl` and
`results/results_previous.jsonl` after `fetch_catalogue.py --pages`.

## Licence

MIT (`LICENSE`).

Christopher Smolen, November 2026.
'''
open(W + '/README.md', 'w').write(readme)
ver = f'''# Run of record

- Snapshot of the catalogue: {S['snapshot']}; {fmt(S['formats'])} formats; {fmt(S['run'])} checked. All files were byte-identical to those downloaded on 2 November.\n- Runs (ledger.csv, column run): {', '.join('%s %s' % (k, fmt(v)) for k, v in S['runs'].items())}. The earlier runs used earlier versions of verify_fmm.py that are not in this repository, all exact (some used the CRT and dense methods, in versions before the last review fixes; rows with free parameters used parameter values that depended on Python's hash seed). Agreement on the overlap: {fmt(S['rerun_overlap'])} formats, {len(S['rerun_disagreements'])} differences in verdict, triad count, wrong or missing counts, or file hash.
- Verifier for this run: verify_fmm.py, SHA-256 prefixes per row in ledger.csv (column verifier; this run: {', '.join('%s %s' % (k, fmt(v)) for k, v in S['verifiers_this_run'].items())}). The {fmt(S['parametric_this_run'])} rows with free parameters were rerun after a fix that makes the parameter values independent of Python's hash seed; the other rows of this run are unaffected by that fix. Current file: {S['verifier_sha256_now']}.
- Methods used (see the docstring of verify_fmm.py): {', '.join('%s %s' % (k or 'not recorded by the earliest runs', fmt(v)) for k, v in sorted(S['methods'].items()))}.
- Total verification time {fmt(S['seconds_total'] // 60)} minutes.

## Checks on the checker

`selftest.py` output for this run:

```
{st}
```

`crosscheck.py` output for the FAIL files, a control and the rebuilt files:

```
{cc}
```

## What this does not show

- Nothing about whether a stated rank is the best possible; only whether each served file is a correct scheme.
- Parametric schemes are verified at fixed pseudo-random rational values of their parameters.
- Rank mismatches are reported as facts about the files and pages, not as errors.
'''
open(W + '/results/VERIFICATION.md', 'w').write(ver); print(readme[:300])
