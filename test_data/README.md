# RAD smoke test

Simulated scTagger long-reads: 98,873 reads, 50 true cells.

Data source: Ebrahimi G, Orabi B, Robinson M, Chauve C, Flannigan R, Hach F.
*Fast and accurate matching of cellular barcodes across short-reads and long-reads
of single-cell RNA-seq experiments.* iScience 25(7):104530 (2022).
DOI: [10.1016/j.isci.2022.104530](https://doi.org/10.1016/j.isci.2022.104530).
Tool: [vpc-ccg/scTagger](https://github.com/vpc-ccg/scTagger).

## Get the data

```bash
curl -L -o test.fq.gz \
  https://github.com/indianewok/rad/releases/download/test-data-v1/shuffled_S_lr_synth.fq.gz
```

True cell barcodes are in each read header: `@<16bp_barcode>-<counter>`.

## scan-wl (de novo cell discovery)

```bash
./rad scan-wl -i test.fq.gz -p CTACACGACGCTCTTCCGATCT -n 16 -w 10x_3v3 -o sw -t 1
```
Expect: `Final barcodes (TRUE): 50`

## demux (correction via 10x_3v3 only)

```bash
./rad demux -l sctagger -q test.fq.gz -d out_default -t 1
```
Expect: ~33,700 reads pass filter, ~71% `CB:Z` exact match to header.

## demux (true_bcs = the 50 ground-truth cells)

```bash
gunzip -c test.fq.gz | awk 'NR%4==1{sub(/^@/,"");sub(/-.*/,"");print}' \
  | sort -u | (echo barcode; cat) > true_50.csv

./rad demux -l sctagger -q test.fq.gz -d out_true -k 10x_3v3 -c true_50.csv -t 1
```
Expect: ~55,700 reads pass filter, **~99.9% `CB:Z` exact match, 50 unique CBs**.
