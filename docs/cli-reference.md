# CLI reference

This is a list of all of the commands and what they get you.

## `rad`

```bash
build/rad <command> [options]
```

Commands:

- `prep` — build/validate a read layout and optional position map
- `scan-wl` — scan reads for the barcodes actually present and build a whitelist
- `demux` — correct barcodes against a whitelist, extract elements, tag reads
- `reformat` — post-process demux output (header rewrite, per-barcode split, coordinates)
- `list` — list registered layout keys and whitelist kit keys with resolved paths
- `modify` — add/remove read-layout or whitelist-kit mappings
- `help` — show help for a specific command

Typical workflow: discover the real barcodes with `scan-wl` (or let `demux -A/--auto-wl`
run it internally), then `demux` against that whitelist, then optional `reformat`. `prep`
is a preliminary layout/position-map step, not part of the per-run barcode path.

Version:

- `rad --version` (also `-V`, `version`) prints the version string.

Exit codes:

- `0` means success
- non-zero means validation/runtime failure

---

## `rad prep`

What it does:

- parses and normalizes layout CSV
- optionally estimates a position map from reads

Usage:

```bash
build/rad prep -l LAYOUT [options]
```

Required:

- `-l, --layout`

Mode flags (set at least one):

- `--read-layout`
- `--position-map`

Options:

- `-q, --fastq` (required with `--position-map`)
- `-o, --output` (required with `--position-map`)
- `-n, --max-reads` (default `50000`)
- `-t, --threads` (default `1`)
- `-v, --verbose`
- `-D, --max-verbose`
- `-h, --help`

Hard checks:

- command fails if neither mode flag is set
- command fails if `--position-map` is set without both `--fastq` and `--output`

Files written:

- with `--position-map`: `<prefix>_layout.csv`, `<prefix>_position_map.csv`
- with `--read-layout` only: output can stay console-only unless `-o` is also set

---

## `rad scan-wl`

What it does:

- scans reads for the sequence sitting immediately downstream of an adapter/primer (the barcode)
- tallies observed barcodes and calls a high-confidence set (the barcodes actually present)
- optionally validates them against a reference whitelist kit/file
- writes a whitelist you then feed to `demux` (or reproduce in one step with `demux -A/--auto-wl`)

This is the barcode-discovery step. Running `demux` against a full reference whitelist
(e.g. the 3M 10x list) corrects reads into any of millions of valid barcodes; `scan-wl`
first narrows that to the cells your data actually contains.

Usage:

```bash
build/rad scan-wl -i FASTQ -p ADAPTER -n BARCODE_LEN -o PREFIX [options]
```

Required:

- `-i, --input` — input FASTQ(.gz) (or `-b, --batch-csv` for batch mode)
- `-p, --adapter_seq` — primer/adapter sequence immediately 5' of the barcode
- `-o, --output-prefix` — output prefix (writes `PREFIX.txt` and `PREFIX.csv`)

Common options:

- `-n, --barcode-length` (default `16`) — bases to extract after the adapter
- `-e, --max-error` (default `0.1`) — adapter max edit-distance ratio
- `-w, --whitelist` — reference kit key or path to validate against (optional)
- `-l, --left-margin` / `-r, --right-margin` (default `0`) — extra bases each side
- `-m, --max-reads` (default all), `-k, --chunk-size` (default `10000`), `-t, --threads`
- `-v, --verbose`, `-h, --help`

Two-part barcode mode (BC1+BC2, e.g. SPLiT-seq / Visium HD): `-1/--bc1-whitelist`,
`-2/--bc2-whitelist`, `-u/--umi-length` (default `9`), `--offset-min`/`--offset-max`.

Files written:

- `<prefix>.txt` — the detected whitelist: high-confidence barcodes, one per line (feed this to `demux -c`)
- `<prefix>.csv` — per-barcode stats (counts + `final_barcode` TRUE/FALSE columns)

---

## `rad demux`

What it does:

- loads layout + position map state
- loads whitelist sets (`true` + `global` model)
- runs chunked `sigalign` extraction/filtering
- writes filtered output and optional debug artifacts

Usage:

```bash
build/rad demux -l LAYOUT -q FASTQ [options]
```

Required:

- `-l, --layout`
- `-q, --fastq`

Whitelist/correction knobs:

- `-k, --kit`
- `-g, --global-whitelist`
- `-c, --custom-whitelist`
- `-R, --bc-correction-mode` (`offensive` default, or `defensive`)
- `--joint-bc-mode` (`default` permissive, or `strict` = every barcode in a joint/split kit must pass)
- `-M, --whitelist-mutation` (default `2`)
- `-m, --generated-mutation` (default `2`)

Auto-whitelist (run `scan-wl` + `demux` in one command):

- `-A, --auto-wl` — before demultiplexing, scan the FASTQ for the barcodes actually present (the `scan-wl` step) and demux against that detected list instead of a full reference. Single-barcode layouts only; split-barcode kits (e.g. Visium HD) should run `rad scan-wl` manually. The adapter and barcode length are derived from the layout; the reference for validation is `-k`/`-g` if given, else the layout's own kit, else de-novo. Writes `<prefix>_scanwl.csv` / `<prefix>_scanwl.txt`.
- `--scan-adapter` — override the derived scan adapter sequence
- `--scan-bc-len` — override the derived barcode length
- `--scan-max-error` (default `0.3`) — adapter max edit-distance ratio for the scan
- `--scan-max-reads` (default `--max-reads`) — reads scanned for the whitelist
- `--scan-chunk` (default `--chunk-size`) — scan chunk size
- `--scan-threads` (default `--threads`) — scan threads

Runtime/output knobs:

- `-n, --max-reads` (default all)
- `-z, --chunk-size` (default `5000`, must be `> 0`)
- `-o, --output` (default `output`)
- `-d, --dir` (default current dir)
- `-F, --log-file`
- `-w, --write-dbg`
- `-b, --bc-split`
- `-t, --threads` (default `1`)
- `--no-umi-rc` — keep UMIs exactly as extracted. By default, for reverse-oriented reads `UB:Z` is reverse-complemented so it shares `CB:Z`'s plus-strand orientation; this flag disables that (raw minus-strand UMIs, e.g. for debugging/metrics).
- `-v, --verbose`
- `-D, --max-verbose`
- `-h, --help`

Practical note:

- `-b/--bc-split` is visible, but split output should be handled by `rad reformat --split-bc` in the current build (the demux-side split is a stub).

Files written:

- primary output: `<outdir>/<prefix>.fq.gz` or `<outdir>/<prefix>.fa.gz` (extension follows input-type logic)
- whitelist summaries: `<outdir>/<prefix>_whitelist_true.csv`, `<outdir>/<prefix>_whitelist_global.csv`
- with `-A/--auto-wl`: `<outdir>/<prefix>_scanwl.csv`, `<outdir>/<prefix>_scanwl.txt` (the internally generated whitelist)
- with `-w`: `<outdir>/<prefix>_dbg.sig.gz`, `<outdir>/<prefix>_dbg.csv.gz`, `<outdir>/<prefix>_dbg.fq.gz` (or `_dbg.fa.gz` — extension follows input-type logic), `<outdir>/<prefix>.metrics.tsv`

Practical note:

- some console path messages may omit `.gz`; those channels are still gzip-compressed on disk.

---

## `rad reformat`

What it does:

- rewrites headers and/or splits reads by barcode tag

Usage:

```bash
build/rad reformat -q INPUT [options]
```

Options:

- `-q, --fastq` (required)
- `-o, --outdir` (required with `--split-bc`)
- `--split-bc` — write one `.fq.gz` per observed `CB:Z` tag
- `--reformat-header` — collapse headers to `QNAME_CB_UB`
- `--reformat-delim CHAR` — delimiter for collapsed headers (single char; `\t`/`\n`/`\r`, or `" "` for space)
- `--coordinate[=MODE]` — rewrite headers using coordinate mode (default MODE `vizHD-v1`; for Visium HD)
- `--bin-size INT` (default `2`) — target bin size in microns for coordinate mode
- `-t, --threads` (default `2`)
- `-v, --verbose`
- `-h, --help`

Hard checks:

- command fails if none of `--split-bc`, `--reformat-header`, or `--coordinate` is set
- command fails if `--split-bc` is set without `--outdir`

Behavior detail:

- `--reformat-header` alone rewrites input in place

Files written:

- `--reformat-header`: input file is rewritten in place
- `--split-bc`: `<outdir>/<CB>.fq.gz` for each observed `CB:Z` tag

---

## `rad list`

Lists registered read-layout keys and whitelist-kit keys with their resolved paths.

```bash
build/rad list                 # both layouts and whitelists
build/rad list --layouts       # layout keys only
build/rad list --whitelists    # whitelist kit keys only
```

Entries that point at a missing file are shown as `[unresolved]` so you can spot broken mappings.

---

## `rad modify`

Adds or removes a read-layout or whitelist-kit mapping.

```bash
build/rad modify --layout KEY --add /path/to/layout.csv     # add/update a layout key
build/rad modify --layout KEY --remove                      # remove (tombstones a built-in)
build/rad modify --whitelist KIT --add /path/to/wl.csv.gz   # add/update a whitelist kit
build/rad modify --whitelist KIT --remove
```

These overrides **persist across invocations**: they are written to
`~/.rad/layout_overrides.tsv` and `~/.rad/whitelist_overrides.tsv` (falling back to
`./.rad_*_overrides.tsv` if `$HOME` is unset). Removing a built-in default is stored as a
tombstone so it stays removed. `build/rad_config` is a standalone legacy helper that writes
the same override files.

---

## Other binaries

These are auxiliary executables built alongside `rad`; the everyday workflow uses the `rad` subcommands above.

### `discover_layout`

```bash
build/discover_layout --help
```

Writes:

- `<prefix>_layout.csv`

### `generate_whitelist`

```bash
build/generate_whitelist --help
```

Standalone build of the `scan-wl` core (prefer `rad scan-wl`). Writes:

- `<output-prefix>.csv`
- `<output-prefix>.txt`

---

## Run checklist

Capture this for every run:

- full command line
- full output file list
- byte size of each major artifact
- read/line counts for main outputs
