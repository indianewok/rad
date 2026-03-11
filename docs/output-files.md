# Output files

## `rad prep`

Typical command:

```bash
build/rad prep -l <layout> --position-map -q <reads> -o <prefix>
```

Writes:

- `<prefix>_layout.csv`
- `<prefix>_position_map.csv`

If only `--read-layout` is used, output may stay console-only unless `-o` is provided.

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

Note:

- some console messages may print paths without `.gz`; actual files are gzip-compressed in those channels.

## `rad reformat`

Header rewrite mode:

```bash
build/rad reformat -q <file.fq.gz> --reformat-header
```

Effect:

- rewrites input in place.

Barcode split mode:

```bash
build/rad reformat -q <file.fq.gz> --split-bc -o <outdir>
```

Writes:

- `<outdir>/<CB>.fq.gz` for each observed `CB:Z` tag.

## `discover_layout`

Writes:

- `<prefix>_layout.csv`

## `generate_whitelist`

Writes:

- `<output-prefix>.csv`
- `<output-prefix>.txt`

## Run checklist

Capture with every run:

- command line used,
- full output file list,
- byte sizes,
- read/line counts for major artifacts.
