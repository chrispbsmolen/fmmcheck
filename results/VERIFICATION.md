# Run of record

- Snapshot of the catalogue: 2026-11-04; 5,222 formats; 5,222 checked. All files were byte-identical to those downloaded on 2 November.
- Runs (ledger.csv, column run): this run 4,101, cloud_2nov 549, cloud_3nov 37, mac_3nov 535. The earlier runs used earlier versions of verify_fmm.py that are not in this repository, all exact (some used the CRT and dense methods, in versions before the last review fixes; rows with free parameters used parameter values that depended on Python's hash seed). Agreement on the overlap: 4,101 formats, 0 differences in verdict, triad count, wrong or missing counts, or file hash.
- Verifier for this run: verify_fmm.py, SHA-256 prefixes per row in ledger.csv (column verifier; this run: b8a06f0e4a7b5f8f 3,581, a561edfa10ca8684 520). The 520 rows with free parameters were rerun after a fix that makes the parameter values independent of Python's hash seed; the other rows of this run are unaffected by that fix. Current file: a561edfa10ca8684.
- Methods used (see the docstring of verify_fmm.py): dense_crt 32, dense_int64 161, none (no scheme) 35, numpy_crt 473, numpy_int64 4,334, python_exact 187.
- Total verification time 228 minutes.

## Checks on the checker

`selftest.py` output for this run:

```
Part A: constructed schemes
naive_2x3x4                        expected PASS               wrong 0    missing 0   ok
strassen                           expected PASS               wrong 0    missing 0   ok
strassen_x_strassen                expected PASS               wrong 0    missing 0   ok
strassen_x_naive_1x2x2             expected PASS               wrong 0    missing 0   ok
strassen_basis_change              expected PASS               wrong 0    missing 0   ok
strassen2_basis_change             expected PASS               wrong 0    missing 0   ok
rank_mismatch                      expected PASS_RANK_MISMATCH wrong 0    missing 0   ok
plus1_first                        expected FAIL               wrong 4    missing 0   ok
plus1_last_triad                   expected FAIL               wrong 16   missing 0   ok
plus_P1                            expected FAIL               wrong 4    missing 0   ok
plus_P1P2                          expected FAIL               wrong 4    missing 0   ok
plus_1_over_P1                     expected FAIL               wrong 16   missing 0   ok
dropped_triad                      expected FAIL               wrong 15   missing 1   ok
empty_triad_added                  expected PASS_RANK_MISMATCH wrong 0    missing 0   ok
huge_denominator                   expected PASS               wrong 0    missing 0   ok
newline_inside_tokens              expected PASS               wrong 0    missing 0   ok
backslash_continuation             expected PASS               wrong 0    missing 0   ok
parametric_l                       expected PASS (parametric)  ok

Part B: catalogue files
catalogue_2x4x5                    expected PASS               wrong 0    missing 0   ok
catalogue_2x4x5_corrupted          expected FAIL               wrong 24   missing 0   ok
2x3x5                              skipped (parametric)
catalogue_2x4x7                    expected PASS               wrong 0    missing 0   ok
catalogue_2x4x7_corrupted          expected FAIL               wrong 6    missing 0   ok
2x3x10                             skipped (parametric)
catalogue_2x3x3                    expected PASS               wrong 0    missing 0   ok
catalogue_2x3x3_corrupted          expected FAIL               wrong 6    missing 0   ok
catalogue_2x3x4                    expected PASS               wrong 0    missing 0   ok
catalogue_2x3x4_corrupted          expected FAIL               wrong 4    missing 0   ok
catalogue_2x3x6                    expected PASS               wrong 0    missing 0   ok
catalogue_2x3x6_corrupted          expected FAIL               wrong 6    missing 0   ok
catalogue_2x4x4                    expected PASS               wrong 0    missing 0   ok
catalogue_2x4x4_corrupted          expected FAIL               wrong 6    missing 0   ok
catalogue_2x2x2                    expected PASS               wrong 0    missing 0   ok
catalogue_2x2x2_corrupted          expected FAIL               wrong 4    missing 0   ok
catalogue_2x3x9                    expected PASS               wrong 0    missing 0   ok
catalogue_2x3x9_corrupted          expected FAIL               wrong 6    missing 0   ok
3x3x3                              skipped (parametric)
catalogue_3x3x6                    expected PASS               wrong 0    missing 0   ok
catalogue_3x3x6_corrupted          expected FAIL               wrong 81   missing 0   ok
catalogue_3x3x5                    expected PASS               wrong 0    missing 0   ok
catalogue_3x3x5_corrupted          expected FAIL               wrong 16   missing 0   ok
catalogue_2x3x7                    expected PASS               wrong 0    missing 0   ok
catalogue_2x3x7_corrupted          expected FAIL               wrong 4    missing 0   ok
catalogue_3x4x5                    expected PASS               wrong 0    missing 0   ok
catalogue_3x4x5_corrupted          expected FAIL               wrong 8    missing 1   ok

ALL TESTS PASS
```

`crosscheck.py` output for the FAIL files, a control and the rebuilt files:

```
data/tensors/8x27x30_tensor.mpl.bz2 triads 3736 trial 0: 96 of 240 entries of A.B wrong mod 2^61-1
data/tensors/15x20x26_tensor.mpl.bz2 triads 4275 trial 0: 366 of 390 entries of A.B wrong mod 2^61-1
data/tensors/26x26x26_tensor.mpl.bz2 triads 8652 trial 0: 676 of 676 entries of A.B wrong mod 2^61-1
data/tensors/28x28x28_tensor.mpl.bz2 triads 10535 trial 0: 784 of 784 entries of A.B wrong mod 2^61-1
data/tensors/30x30x30_tensor.mpl.bz2 triads 12672 trial 0: 900 of 900 entries of A.B wrong mod 2^61-1
data/tensors/32x32x32_tensor.mpl.bz2 triads 15079 trial 0: 1024 of 1024 entries of A.B wrong mod 2^61-1
data/tensors/12x28x28_tensor.mpl.bz2 triads 5328 trial 0: 0 of 336 entries of A.B wrong mod 2^61-1
data/fixes/11x25x28_tensor_rebuilt.mpl.bz2 triads 4572 trial 0: 0 of 308 entries of A.B wrong mod 2^61-1
data/fixes/15x20x26_tensor_rebuilt.mpl.bz2 triads 4466 trial 0: 0 of 390 entries of A.B wrong mod 2^61-1
data/fixes/8x27x30_tensor_rebuilt.mpl.bz2 triads 3750 trial 0: 0 of 240 entries of A.B wrong mod 2^61-1
```

## What this does not show

- Nothing about whether a stated rank is the best possible; only whether each served file is a correct scheme.
- Parametric schemes are verified at fixed pseudo-random rational values of their parameters.
- Rank mismatches are reported as facts about the files and pages, not as errors.
