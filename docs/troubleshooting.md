# Troubleshooting

Format here is simple: symptom -> why -> what to do.

## OpenMP build failure

Symptom:
- OpenMP not found or linker OpenMP errors.

Why:
- compiler/runtime mismatch (especially macOS default clang setups).

What to do:
- install `llvm` + `libomp`.
- configure CMake with LLVM compilers.
- rebuild cleanly.

See [`installation.md`](installation.md).

## `Unknown layout key`

Symptom:
- `prep`/`demux` fails on unknown layout.

Why:
- key not present in runtime layout map.

What to do:

```bash
build/rad_config layout list
```

Use one of those keys, or pass an absolute layout CSV path.

## `Cannot locate 'resources' folder`

Symptom:
- runtime exception for missing `resources`.

Why:
- executable moved away from expected relative resource path.

What to do:
- run from repository tree, or
- deploy `resources/` with the binary path layout RAD expects.

## `Whitelist file not found`

Symptom:
- whitelist load failure.

Why:
- alias/path resolution failed (bad key/path/whitespace).

What to do:

```bash
build/rad_config whitelist list
```

Then use explicit absolute paths if needed.

## Slow `.gz` throughput

Symptom:
- compressed I/O is much slower than expected.

Why:
- fallback gzip path rather than `pigz` path.

What to do:
- install `pigz`.
- set `RAD_PIGZ` if needed.
- tune `RAD_PIGZ_THREADS`.

## Empty or near-empty demux output

Symptom:
- output FASTQ exists but very few reads pass.

Why (usual causes):
- layout/chemistry mismatch,
- whitelist mismatch,
- overly strict correction assumptions.

What to do:
- inspect layout with `prep --read-layout`.
- rerun with `-v` and `-w`.
- inspect `*_dbg.csv.gz` + `*.metrics.tsv`.
- verify kit/whitelist pairing.

## Output order changes b/w runs

Symptom:
- same input, different record order across runs.

Why:
- multi-threaded chunk scheduling and async write behavior.

What to do:
- use `--threads 1` when stable ordering is required.
- compare canonicalized content rather than raw line order for multi-thread runs.

## Behavior looks stale after edits

Symptom:
- runtime/help output does not match recent code changes.

Why:
- stale build artifacts.

What to do:

```bash
cmake --build build --clean-first
```
