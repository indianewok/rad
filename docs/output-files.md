# Output files

This page is the file-level contract: what each command writes and where it lands.

## `rad prep`

Typical command:

```bash
build/rad prep -l <layout> --position-map -q <reads> -o <prefix>
```

Writes:
- `<prefix>_layout.csv`
- `<prefix>_position_map.csv`

If you only run `--read-layout`, output can stay console-only unless you also pass `-o`.

## `rad demux`

Typical command:

```bash
build/rad demux -l <layout> -q <reads> -o <prefix> -d <outdir>
```

Primary output:
- `<outdir>/<prefix>.fq.gz` or `<outdir>/<prefix>.fa.gz` (extension follows input type logic)

Whitelist summaries:
- `<outdir>/<prefix>_whitelist_true.csv`
- `<outdir>/<prefix>_whitelist_global.csv`

With `-w` debug enabled:
- `<outdir>/<prefix>_dbg.sig.gz`
- `<outdir>/<prefix>_dbg.csv.gz`
- `<outdir>/<prefix>_dbg.fq.gz`
- `<outdir>/<prefix>.metrics.tsv`

Practical note:
- some console path messages may omit `.gz`; those channels are still gzip-compressed on disk.

## `rad reformat`

Header rewrite mode:

```bash
build/rad reformat -q <file.fq.gz> --reformat-header
```

Effect:
- rewrites the input file in place

Barcode split mode:

```bash
build/rad reformat -q <file.fq.gz> --split-bc -o <outdir>
```

Writes:
- `<outdir>/<CB>.fq.gz` for each observed `CB:Z` tag

## `discover_layout`

Writes:
- `<prefix>_layout.csv`

## `generate_whitelist`

Writes:
- `<output-prefix>.csv`
- `<output-prefix>.txt`

## Run checklist

Capture this with every run:
- full command line
- full output file list
- byte size of each major artifact
- read/line counts for the main outputs
