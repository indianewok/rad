# RAD docs index

## Core usage

- [`installation.md`](installation.md): build from source or install via Bioconda (`mamba install -c bioconda rad`), plus the first run — discover the real barcodes with `scan-wl` (or `demux -A/--auto-wl`), then `demux`, then optional `reformat`.
- [`cli-reference.md`](cli-reference.md): every command (`prep`, `scan-wl`, `demux`, `reformat`, `list`, `modify`), its flag contract, and exact output files.
- [`layouts-and-whitelists.md`](layouts-and-whitelists.md): layout CSV format, bundled layouts, whitelist sizes/sources, kit pairings, and building a whitelist from your reads with `scan-wl`.

## Internals

- [`architecture.md`](architecture.md): how layout parsing, whitelist correction, and chunked processing fit together.
