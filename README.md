# fmmcheck

An exact checker for Alexandre Sedoglavic's catalogue of fast matrix multiplication algorithms at
https://fmm.univ-lille.fr/, and the results of running it on every scheme file in the catalogue.

Each catalogue entry `<X×Y×Z:r>` has a file `XxYxZ_tensor.mpl.bz2` holding triads `(U_i, V_i, W_i)` and the
identity `A.B = sum_i Tr(U_i^T A) Tr(V_i^T B) W_i^T`. `verify_fmm.py` proves or refutes that identity in exact
arithmetic (for schemes with free parameters, at fixed parameter values), by evaluating every coefficient of the
Brent equations, and reports whether the number of triads equals the rank in the catalogue's table. Two cheaper
checks multiply random matrices with the file's own formula: `randeval.py` (exact rationals, sharing only the file
parser with `verify_fmm.py`) and `crosscheck.py` (written separately, with its own parser, modulo 2^61 - 1).

## Results, snapshot of 2026-11-04

All 5,222 formats in the front-page table were checked (0 not run). 4,101 were checked in this run with the code in this repository. The other 1,121, mostly the largest formats and including the four square formats below, were checked only in earlier runs of 2 and 3 November with earlier versions of `verify_fmm.py` that are not in this repository (also exact; the ledger's run column says which). On the 4,101 formats checked in both, the results agree exactly (0 differences).

| verdict | formats |
|---|---|
| PASS: the equations hold and the triad count equals the table's rank | 4,001 |
| PASS_RANK_MISMATCH: the equations hold, but the file has a different number of triads than the table | 1,180 |
| FAIL: the equations do not hold | 6 |
| no scheme in the file | 35 |

Rank mismatches, compared with the rank stated on each format's own page (`results/rank_mismatches.csv`):
855 where the format page agrees with the file and only the front page differs,
195 where the format page agrees with the front page and the file differs,
87 where the format page is empty or gives no rank and 43 where it differs from both.
In 11 of them the file has fewer triads than the table states.

Files with no scheme: 30 placeholder (Tensor:=Tensor:, no triads), 4 stub (calls ReturnFastMatrixMultiplication, no triads), 1 truncated file (text ends inside a matrix).

### FAIL

| format | front page | format page | file triads | wrong / missing Brent coefficients | front-page source |
|---|---|---|---|---|---|
| 8x27x30 | 3744 | 3744 | 3736 | 28062 wrong, 36 missing | ( ⟨4×9×10:250⟩ - 12) ⊗ ⟨2×3×3:15⟩ +6 ⟨4×3×3:29⟩ |
| 15x20x26 | 4464 | 4464 | 4275 | 111480 wrong, 936 missing | ( ⟨5×5×13:227⟩ - 64) ⊗ ⟨3×4×2:20⟩ +4 ⟨3×4×6:54⟩ +26 ⟨3×4×4:38⟩ |
| 26x26x26 | 8051 | 8652 | 8652 | 302673735 wrong, 0 missing | Khoruzhii, Gelß and Pokutta 2026 [  GitHub repository  ] |
| 28x28x28 | 9846 | 10535 | 10535 | 472423014 wrong, 0 missing | Khoruzhii, Gelß and Pokutta 2026 [  GitHub repository  ] |
| 30x30x30 | 11889 | 12672 | 12672 | 715016295 wrong, 0 missing | Khoruzhii, Gelß and Pokutta 2026 [  GitHub repository  ] |
| 32x32x32 | 14196 | 15079 | 15079 | 1053562608 wrong, 0 missing | Khoruzhii, Gelß and Pokutta 2026 [  GitHub repository  ] |

Every FAIL is confirmed by `crosscheck.py` (`results/crosscheck_output.txt`): random products are wrong in many
entries, while a correct file of similar size passes. For the four square formats the front page credits Khoruzhii,
Gelß and Pokutta 2026 at the ranks above and links https://github.com/khoruzhii/lita , while the served files match
the format pages (rank 8652 and so on, "private communication june 2026"). The schemes in that repository (commit
472a79f, `git clone https://github.com/khoruzhii/lita`), at the front-page ranks, were checked here:

| lita scheme (commit 472a79f) | rank | exact random-matrix check (tools/lita_randeval.py) | exact proof (tools/lita_exact.py) |
|---|---|---|---|
| 26x26x26_r8051.npz | 8051 | 0 of 676 entries wrong | started, not completed |
| 28x28x28_r9846.npz | 9846 | 0 of 784 entries wrong | not run (too slow on the machine used) |
| 30x30x30_r11889.npz | 11889 | 0 of 900 entries wrong | not run (too slow on the machine used) |
| 32x32x32_r14196.npz | 14196 | 0 of 1024 entries wrong | not run (too slow on the machine used) |

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

## Credit, corrections and removal

Every result here names its source; please cite the original authors for their results. Corrections, missing or
wrong credit, and requests to change or remove anything (including a request to stop linking a file) are welcome as
a GitHub issue or by email to the maintainer, Christopher Smolen (independent researcher, ORCID 0009-0002-9537-7961),
chrispbsmolen@gmail.com. A file whose author objects to its inclusion is removed. Every such change is recorded,
with its date and reason, in `CORRECTIONS.md`.

## Licence

MIT (`LICENSE`).

Christopher Smolen, November 2026.
