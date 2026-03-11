# CLI reference

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
- `0` success
- non-zero validation/runtime failure

---

## `rad prep`

What it does:
- parses/normalizes layout,
- optionally estimates position map from reads.

Usage:

```bash
build/rad prep -l LAYOUT [options]
```

Required:
- `-l, --layout`

Mode flags (at least one required):
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
- fails if neither mode flag is set.
- fails if `--position-map` is set without `--fastq` and `--output`.

---

## `rad demux`

What it does:
- loads layout + whitelist state,
- runs `sigalign` pipeline,
- writes filtered output plus optional debug artifacts.

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
- `-b, --bc_split` (see caveat)
- `-t, --threads` (default `1`)
- `-v, --verbose`
- `-D, --max_verbose`
- `-h, --help`

Current caveat:
- `--bc_split` is exposed in help, but split execution inside `demux` is currently stubbed.
- Use `rad reformat --split-bc` for actual split output.

---

## `rad reformat`

What it does:
- rewrites headers and/or splits reads by barcode tag.

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
- fails if neither `--split-bc` nor `--reformat-header` is set.
- fails if `--split-bc` is set without `--outdir`.

Behavior detail:
- `--reformat-header` alone rewrites input in place.

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

Current caveat:
- `set/rm` updates are process-local and not persistent across independent invocations.

---

## Other binaries

### `discover_layout`

```bash
build/discover_layout --help
```

### `generate_whitelist`

```bash
build/generate_whitelist --help
```
