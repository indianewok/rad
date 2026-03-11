# Troubleshooting

Format is simple: symptom -> why it happens -> what to do next.

## OpenMP build failure

Symptom:
- OpenMP isn't found, or linking fails with OpenMP errors

Why:
- compiler/runtime mismatch (especially default macOS clang setups)

What to do:
- install `llvm` + `libomp`
- configure CMake with LLVM compilers
- rebuild cleanly

See [`installation.md`](installation.md).

## `Unknown layout key`

Symptom:
- `prep` or `demux` fails on an unknown layout

Why:
- the key isn't in the runtime layout map

What to do:

```bash
build/rad_config layout list
```

Use one of those keys, or pass an absolute layout CSV path.

## `Cannot locate 'resources' folder`

Symptom:
- runtime exception says `resources` can't be found

Why:
- binary moved away from the expected relative resource path

What to do:
- run from the repository tree, or
- deploy `resources/` with the binary path layout RAD expects

## `Whitelist file not found`

Symptom:
- whitelist load fails

Why:
- alias/path resolution failed (bad key/path/whitespace)

What to do:

```bash
build/rad_config whitelist list
```

Then use explicit absolute paths if needed.

## Slow `.gz` throughput

Symptom:
- compressed I/O is much slower than expected

Why:
- RAD fell back to the non-`pigz` path

What to do:
- install `pigz`
- set `RAD_PIGZ` if needed
- tune `RAD_PIGZ_THREADS`

## Empty or near-empty demux output

Symptom:
- output FASTQ exists, but almost no reads pass

Why (common causes):
- layout/chemistry mismatch
- whitelist mismatch
- correction settings are too strict

What to do:
- inspect layout with `prep --read-layout`
- rerun with `-v` and `-w`
- inspect `*_dbg.csv.gz` and `*.metrics.tsv`
- verify kit/whitelist pairing

## Output order changes between runs

Symptom:
- same input, different record order across runs

Why:
- chunk scheduling and async writes in multithreaded mode

What to do:
- use `--threads 1` if stable ordering matters
- compare canonicalized content, not raw line order, for multithreaded runs

## Behavior looks stale after edits

Symptom:
- runtime/help output doesn't match recent code changes

Why:
- stale build artifacts

What to do:

```bash
cmake --build build --clean-first
```
