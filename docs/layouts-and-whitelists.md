# Layouts and whitelists

This page is the practical map: what’s bundled, what it maps to, and what to try first.

## Bundled layouts in the current build

```bash
build/rad_config layout list
```

| Layout key | Source file | Typical pattern | Default whitelist field in layout |
| --- | --- | --- | --- |
| `five_prime` | `resources/read_layout/five_prime_read_layout.csv` | 10x-like 5' (primer, barcode, UMI, TSO, read) | `10x_5v1` |
| `three_prime` | `resources/read_layout/three_prime_read_layout.csv` | 10x-like 3' (primer, barcode, UMI, polyT, read) | `10x_3v3` |
| `sctagger` | `resources/read_layout/sctagger_sim_read_layout.csv` | simulated/tagging 3'-style layout | `10x_3v3` |
| `splitseq` | `resources/read_layout/splitseq_read_layout.csv` | SPLiT-seq multi-round barcodes + linkers | `splitseq_bc1`, `splitseq_bc2` |
| `visium` | `resources/read_layout/visium_three_prime_read_layout.csv` | Visium-like spatial 3' layout | `10x_Vis_V1` |
| `nanopore_rapid_bc` | `resources/read_layout/nanopore_bulk_rapid_bc_read_layout.csv` | Nanopore bulk rapid barcode-like template | none |

Where these come from:
- all are local templates shipped in `resources/read_layout/`.
- nothing is downloaded at runtime.

Visium HD status note:
- Visium HD support is being developed on `dev`.
- In the current `main` build/docs format, treat Visium HD as not fully integrated into the stable layout/whitelist flow yet.

## Quick layout pick guide

- 10x 3' style library -> start with `three_prime`
- 10x 5' style library -> start with `five_prime`
- SPLiT-seq -> start with `splitseq`
- Visium -> start with `visium`
- unclear chemistry -> pick closest, run `prep --read-layout`, then edit a custom CSV

## Custom layout CSV schema

Columns RAD expects:
- `id`
- `seq`
- `expected_length`
- `type`
- `class`
- `whitelist`
- optional `direction`

Minimal example:

| id | seq | expected_length | type | class | whitelist |
| --- | --- | --- | --- | --- | --- |
| `forw_primer` | `CTACACGACGCTCTTCCGATCT` |  | `static` |  |  |
| `barcode` |  | `16` | `variable` | `barcode` | `10x_5v1` |
| `umi` |  | `10` | `variable` | `umi` |  |
| `read` |  |  | `variable` | `read` |  |
| `rev_primer` | `GTACTCTGCGTTGATACCACTGCTT` |  | `static` |  |  |

Direction notes:
- if `direction` is omitted, RAD treats rows as forward and auto-generates reverse-complement counterparts.
- `forward_only` / `reverse_only` can be used for one-sided elements.

## Direction column examples (`both`, `forward_only`, `reverse_only`)

### 1) Both orientations (auto forward + reverse-complement)

Use the standard header and omit `direction`:

```csv
id,seq,expected_length,type,class,whitelist
forw_primer,CTACACGACGCTCTTCCGATCT,,static,,
barcode,,16,variable,barcode,10x_5v1
```

### 2) Forward-only element(s)

Add a `direction` column and mark specific rows as `forward_only`:

```csv
id,seq,expected_length,type,class,whitelist,direction
barcode,,16,variable,barcode,10x_5v1,forward_only
read,,,variable,read,,forward_only
```

### 3) Reverse-only element(s)

Use the same header and mark rows as `reverse_only`:

```csv
id,seq,expected_length,type,class,whitelist,direction
adapter_rc,AGATCGGAAGAGCGTCGTGTAG,,static,adapter,,reverse_only
barcode_rc,,16,variable,barcode,10x_5v1,reverse_only
```

## Bundled whitelist resources (current files + sizes)

```bash
build/rad_config whitelist list
```

Entry counts below are line counts from bundled `.gz` files.

| File | Entries | Typical source/meaning |
| --- | ---: | --- |
| `737K-august-2016_bitlist.csv.gz` | 737,280 | 10x 737K-era barcode family |
| `3M-february-2018-3v3.txt_bitlist.csv.gz` | 6,794,880 | 10x 3' v3/v3.1 family |
| `3M-3pgex-may-2023.txt_bitlist.csv.gz` | 7,372,800 | 10x 3' newer chemistry family |
| `3M-5pgex-jan-2023.txt_bitlist.csv.gz` | 3,686,400 | 10x 5' v3 family |
| `visium-v1_v2_bitlist.csv.gz` | 4,992 | Visium v1/v2 spatial barcodes |
| `visium-v3_v4_bitlist.csv.gz` | 4,993 | Visium v3/v4 spatial barcodes |
| `visium-v5_bitlist.csv.gz` | 14,337 | Visium v5 spatial barcodes |
| `splitseq_bc1_bitlist.csv.gz` | 96 | SPLiT-seq barcode set |
| `splitseq_bc2_bitlist.csv.gz` | 96 | SPLiT-seq barcode set |
| `visium_hd_coordinates.csv.gz` | 11,222,501 | Visium HD coordinate table (auxiliary spatial resource) |

## How RAD stores barcodes internally (`int64_seq` + `barcode_entry`)

RAD stores barcodes in a compact 2-bit encoding (`int64_seq`), then wraps that in `barcode_entry` for counts/flags.

### `int64_seq` structure

- `length` (`uint16_t`): barcode length in bases.
- `bits` (`vector<int64_t>`): packed 2-bit chunks (up to 32 bases per `int64_t` chunk).

Base map used in encoding:

| Base | 2-bit code |
| --- | --- |
| `A` | `00` |
| `C` | `01` |
| `T` | `10` |
| `G` | `11` |

Encoding rule per chunk:
- initialize `result = 0`
- for each base (left to right): `result = (result << 2) | code(base)`

That means:
- first base in a chunk is stored in the highest used bits of that chunk
- last base in a chunk is stored in the least-significant 2 bits

Important positional detail:
- many low-level mutation/distance routines index from LSB upward (`(bits >> (2*pos)) & 3`),
- so `pos = 0` corresponds to the **rightmost** base in the decoded sequence.

Character constraints:
- encoder accepts `A/C/T/G` only (uppercase in the direct `int64_seq` conversion path).
- `N` is not encodable in this 2-bit schema.

### `barcode_entry` structure

- `barcode` (`int64_seq`)
- `count` (`counter`, atomic count array)
- `filtered` (`bool`)

Counter categories are:
- `raw`, `forw`, `forw_concat`, `rev`, `rev_concat`, `total`, `corrected`, `filtered`

### Translating bitlist values to DNA barcodes

Bitlist whitelist files store integer-encoded barcodes.  
To decode correctly, the barcode length is required (the integer alone is not enough).

Python decode helper:

```python
def decode_barcode(value: int, length: int) -> str:
    # if negative values ever appear, reinterpret as unsigned 64-bit
    if value < 0:
        value &= (1 << 64) - 1

    lut = "ACTG"  # 0->A, 1->C, 2->T, 3->G
    out = []
    for _ in range(length):
        out.append(lut[value & 0b11])
        value >>= 2
    return "".join(reversed(out))
```

Python encode helper:

```python
def encode_barcode(seq: str) -> int:
    code = {"A": 0, "C": 1, "T": 2, "G": 3}
    v = 0
    for b in seq.upper():
        v = (v << 2) | code[b]
    return v
```

Example:
- `decode_barcode(23904338, 16)` -> `AAACCTGAGAAACCAT`

## Kit aliases -> whitelist files (current map)

| Kit group | Aliases | Backing file |
| --- | --- | --- |
| 10x 3' older | `10x_3v1`, `10x_3v2` | `737K-august-2016_bitlist.csv.gz` |
| 10x 3' v3 family | `10x_3v3`, `10x_3v3.1`, `10x_3HTv3.1` | `3M-february-2018-3v3.txt_bitlist.csv.gz` |
| 10x 3' newer | `10x_3v4` | `3M-3pgex-may-2023.txt_bitlist.csv.gz` |
| 10x 5' older | `10x_5v1`, `10x_5v2`, `10x_5HTv2` | `737K-august-2016_bitlist.csv.gz` |
| 10x 5' newer | `10x_5v3` | `3M-5pgex-jan-2023.txt_bitlist.csv.gz` |
| Visium v1/v2 | `10x_Vis_V1`, `10x_Vis_V2` | `visium-v1_v2_bitlist.csv.gz` |
| Visium v3/v4 | `10x_Vis_V3`, `10x_Vis_V4` | `visium-v3_v4_bitlist.csv.gz` |
| Visium v5 | `10x_Vis_V5` | `visium-v5_bitlist.csv.gz` |
| SPLiT-seq | `splitseq_bc1`, `splitseq_bc2` | `splitseq_bc1_bitlist.csv.gz`, `splitseq_bc2_bitlist.csv.gz` |

## Common pairings (quick start)

- `three_prime` + `10x_3v3` (or the specific 3' alias for the kit)
- `five_prime` + `10x_5v1`/`10x_5v2`/`10x_5v3` (chemistry-dependent)
- `visium` + `10x_Vis_V1`..`10x_Vis_V5`
- `splitseq` + `splitseq_bc1` + `splitseq_bc2`

If the best pairing is unclear, run a short pilot with `-w` and inspect debug outputs before scaling up.

## Custom whitelist files

Explicit file paths always work and are usually the safest in production scripts:

```bash
build/rad demux \
  -l five_prime \
  -q reads.fq.gz \
  -g /data/global_wl.csv.gz \
  -c /data/true_wl.csv.gz \
  -o demo -d run
```

## `rad_config` caveat

`rad_config set/rm` is process-local in the current implementation, so it does not persist across separate invocations.
