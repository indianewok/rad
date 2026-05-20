# CLI reference

This is a list of all of the commands and what they get you.

## `rad`

```bash
build/rad <command> [options]
```

Commands:

- `prep`
- `demux`
- `reformat`
- `help`

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
- `-n, --max_reads` (default `50000`)
- `-t, --threads` (default `1`)
- `-v, --verbose`
- `-D, --max_verbose`
- `-h, --help`

Hard checks:

- command fails if neither mode flag is set
- command fails if `--position-map` is set without both `--fastq` and `--output`

Files written:

- with `--position-map`: `<prefix>_layout.csv`, `<prefix>_position_map.csv`
- with `--read-layout` only: output can stay console-only unless `-o` is also set

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
- `-g, --global_whitelist`
- `-c, --custom_whitelist`
- `-R, --bc_correction_mode` (`offensive` default, or `defensive`)
- `-M, --whitelist_mutation` (default `2`)
- `-m, --generated_mutation` (default `2`)

Runtime/output knobs:

- `-n, --max_reads` (default all)
- `-z, --chunk_size` (default `5000`, must be `> 0`)
- `-o, --output` (default `output`)
- `-d, --dir` (default current dir)
- `-F, --log-file`
- `-w, --write_dbg`
- `-b, --bc_split`
- `-t, --threads` (default `1`)
- `-v, --verbose`
- `-D, --max_verbose`
- `-h, --help`

Practical note:

- `--bc_split` is visible, but split output should be handled by `rad reformat --split-bc` in the current build.

Files written:

- primary output: `<outdir>/<prefix>.fq.gz` or `<outdir>/<prefix>.fa.gz` (extension follows input-type logic)
- whitelist summaries: `<outdir>/<prefix>_whitelist_true.csv`, `<outdir>/<prefix>_whitelist_global.csv`
- with `-w`: `<outdir>/<prefix>_dbg.sig.gz`, `<outdir>/<prefix>_dbg.csv.gz`, `<outdir>/<prefix>_dbg.fq.gz`, `<outdir>/<prefix>.metrics.tsv`

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
- `--split-bc`
- `--reformat-header`
- `-t, --threads` (default `2`)
- `-v, --verbose`
- `-h, --help`

Hard checks:

- command fails if neither `--split-bc` nor `--reformat-header` is set
- command fails if `--split-bc` is set without `--outdir`

Behavior detail:

- `--reformat-header` alone rewrites input in place

Files written:

- `--reformat-header`: input file is rewritten in place
- `--split-bc`: `<outdir>/<CB>.fq.gz` for each observed `CB:Z` tag

---

## `rad_config`

Usage:

```bash
build/rad_config layout list
build/rad_config layout get <type>
build/rad_config layout set <type> <path>
build/rad_config layout rm <type>

build/rad_config whitelist list
build/rad_config whitelist get <kit>
build/rad_config whitelist set <kit> <path>
build/rad_config whitelist rm <kit>
```

Practical note:

- `set/rm` updates are process-local in the current build and won't persist across independent invocations.

---

## Other binaries

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

Writes:

- `<output-prefix>.csv`
- `<output-prefix>.txt`

---

## Run checklist

Capture this for every run:

- full command line
- full output file list
- byte size of each major artifact
- read/line counts for main outputs
