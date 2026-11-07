#!/usr/bin/env python3
"""Exercise the demux stdout interface with a self-contained cached layout."""
import csv
import gzip
import os
from pathlib import Path
import random
import select
import signal
import subprocess
import sys
import tempfile


FORWARD = "CTTCCGATCTATGGCGACCTTATCAG"
LINKER = "TTGTCTTCCTAAGAC"
REVERSE = "CCCGCCTCTCAGTACGTCAGCAG"
BARCODES = ["ACGTCAGTACGATCGTAGCTACGTA", "GCTAGCATCGTACAGTCGATGCTAC"]
LAYOUT_HEADER = ["id", "seq", "masked_seq", "expected_length", "length_candidates", "flags",
                 "type", "class", "direction", "class_id", "whitelist", "order"]
MAP_HEADER = ["id", "primary_start", "primary_stop", "secondary_start", "secondary_stop",
              "misalign_lower", "misalign_mean", "misalign_upper", "align_start", "var_start",
              "align_stop", "var_stop", "misalign_start", "mvar_start", "misalign_stop", "mvar_stop"]


def require(ok, message):
    if not ok:
        raise AssertionError(message)


def revcomp(text):
    return text.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def fixture(root):
    whitelist = root / "barcodes.txt"
    whitelist.write_text("\n".join(BARCODES) + "\n")
    rows = [
        ("seq_start", "", 0, "static", "start", "forward"),
        ("forw_primer", FORWARD, len(FORWARD), "static", "forw_primer", "forward"),
        ("barcode", "", 25, "variable", "barcode", "forward"),
        ("linker", LINKER, len(LINKER), "static", "linker", "forward"),
        ("umi", "", 10, "variable", "umi", "forward"),
        ("poly_t", "T{12,}+", 12, "static", "poly_tail", "forward"),
        ("read", "", "", "variable", "read", "forward"),
        ("rev_primer", REVERSE, len(REVERSE), "static", "rev_primer", "forward"),
        ("seq_stop", "", 0, "static", "stop", "forward"),
        ("rc_seq_start", "", 0, "static", "start", "reverse"),
        ("rc_rev_primer", revcomp(REVERSE), len(REVERSE), "static", "rev_primer", "reverse"),
        ("rc_read", "", "", "variable", "read", "reverse"),
        ("poly_a", "A{12,}+", 12, "static", "poly_tail", "reverse"),
        ("rc_umi", "", 10, "variable", "umi", "reverse"),
        ("rc_linker", revcomp(LINKER), len(LINKER), "static", "linker", "reverse"),
        ("rc_barcode", "", 25, "variable", "barcode", "reverse"),
        ("rc_forw_primer", revcomp(FORWARD), len(FORWARD), "static", "forw_primer", "reverse"),
        ("rc_seq_stop", "", 0, "static", "stop", "reverse"),
    ]
    layout = root / "fixture_layout.csv"
    with layout.open("w", newline="") as handle:
        handle.write("Read Layout,,,,,,,,,,,\n")
        writer = csv.writer(handle)
        writer.writerow(LAYOUT_HEADER)
        for order, (name, seq, length, typ, klass, direction) in enumerate(rows, 1):
            masked = seq if typ == "static" and klass not in ("start", "stop", "poly_tail") else ""
            writer.writerow([name, seq, masked, length, length if length else "", "", typ,
                             klass, direction, name, str(whitelist) if klass == "barcode" else "", order])
    refs = {
        "barcode": ["forw_primer|stop+1", "forw_primer|stop+25", "linker|start-25", "linker|start-1"],
        "umi": ["linker|stop+1", "linker|stop+10", "linker|stop+1|left_terminal_linked", "linker|stop+10|left_terminal_linked"],
        "read": ["poly_t|stop+1", "rev_primer|start-1", "umi|stop+1|left_truncated", "seq_stop|start-1|right_truncated"],
        "rc_read": ["rc_rev_primer|stop+1", "poly_a|start-1", "rc_seq_start|stop+1|left_truncated", "rc_umi|start-1|right_truncated"],
        "rc_umi": ["rc_linker|start-10", "rc_linker|start-1", "rc_linker|start-10|right_terminal_linked", "rc_linker|start-1|right_terminal_linked"],
        "rc_barcode": ["rc_forw_primer|start-25", "rc_forw_primer|start-1", "rc_linker|stop+1", "rc_linker|stop+25"],
    }
    with (root / "fixture_position_map.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(MAP_HEADER)
        for name, seq, _, typ, klass, _ in rows:
            # The fixture has perfect adapters, so do not let the fallback
            # permissive threshold find extra adapter fragments in its payload.
            # Broad cached position bounds keep both read orientations valid.
            stats = ([0, 0, 0, 50, 100, 100, 100, 50, 100, 100, 100]
                     if typ == "static" and seq and klass != "poly_tail" else [""] * 11)
            writer.writerow([name, *refs.get(name, [""] * 4), *stats])
    generator = random.Random(183431)
    units = []
    for i in range(6):
        payload = "".join(generator.choice("ACGT") for _ in range(90 + i * 7))
        umi = "".join(generator.choice("ACGT") for _ in range(10))
        sequence = FORWARD + BARCODES[i % 2] + LINKER + umi + "T" * 16 + payload + REVERSE
        if i % 3 == 1:
            sequence = revcomp(sequence)
        units.append(sequence)
    fastq = root / "input.fastq"
    fasta = root / "input.fasta"
    write_reads(fastq, units)
    write_reads(fasta, units, fasta=True)
    return layout, fastq, fasta, units


def write_reads(path, sequences, fasta=False):
    with path.open("w") as handle:
        for i, sequence in enumerate(sequences):
            if fasta:
                handle.write(">read{} original-comment\n{}\n".format(i, sequence))
            else:
                quality = "".join(chr(33 + j % 50) for j in range(len(sequence)))
                handle.write("@read{} original-comment\n{}\n+\n{}\n".format(i, sequence, quality))


def parse(data, fasta=False):
    require(not data.startswith(b"\x1f\x8b"), "stdout was gzip-compressed")
    lines = data.splitlines()
    step = 2 if fasta else 4
    require(len(lines) % step == 0, "stream contains an incomplete sequence record or diagnostic bytes")
    records = []
    for start in range(0, len(lines), step):
        row = lines[start:start + step]
        require(row[0].startswith(b">" if fasta else b"@"), "diagnostic text contaminated sequence stdout")
        require(b"CB:Z:" in row[0] and b"UB:Z:" in row[0], "barcode or UMI tag missing from stdout")
        require(row[1] and set(row[1]) <= set(b"ACGTN"), "invalid sequence in stdout")
        if not fasta:
            require(row[2].startswith(b"+") and len(row[1]) == len(row[3]), "FASTQ sequence/quality mismatch")
        records.append(tuple(row))
    return records


def no_stream_artifacts(folder):
    forbidden = []
    for path in folder.rglob("*"):
        if not path.is_file():
            continue
        name = path.name
        if "_dbg" in name or name.endswith(".metrics.tsv") or any(
                name.endswith(suffix) for suffix in (".fq", ".fq.gz", ".fastq", ".fastq.gz", ".fa", ".fa.gz", ".fasta", ".fasta.gz")):
            forbidden.append(str(path))
    require(not forbidden, "stdout mode wrote sequence/debug artifacts: " + repr(forbidden))


def main():
    require(len(sys.argv) == 2, "usage: stdout_cli_test.py /path/to/rad")
    binary = str(Path(sys.argv[1]).resolve())
    with tempfile.TemporaryDirectory(prefix="rad-stdout-cli-") as temporary:
        root = Path(temporary)
        layout, fastq, fasta, units = fixture(root)
        env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")

        def command(folder, source=fastq, extra=()):
            folder.mkdir()
            return [binary, "demux", "--layout", str(layout), "--fastq", str(source),
                    "--dir", str(folder), "--output", "result", "--threads", "1",
                    "--chunk-size", "2", "--min-read-length", "0", *extra]

        def run(name, extra=(), source=fastq, expected=0, local_env=None):
            folder = root / name
            args = command(folder, source, extra)
            result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    env=local_env or env, timeout=35)
            require(result.returncode == expected,
                    "{} returned {} (expected {})\n{}".format(name, result.returncode, expected,
                                                               result.stderr.decode(errors="replace")[-12000:]))
            return folder, result

        folder, normal = run("normal")
        default_path = folder / "result.fq.gz"
        require(default_path.exists(), "default mode stopped writing its FASTQ file")
        baseline = gzip.decompress(default_path.read_bytes())
        expected_records = parse(baseline)
        require(len(expected_records) == len(units), "synthetic fixture did not produce all intended reads")
        require(b"[main]" in normal.stdout, "default diagnostics changed unexpectedly")

        folder, streamed = run("stdout", ["--stdout"])
        require(parse(streamed.stdout) == expected_records, "stdout records differ from normal output")
        require(b"[main]" in streamed.stderr, "streaming diagnostics did not reach stderr")
        no_stream_artifacts(folder)

        for name, extra in [("verbose", ["-v"]), ("max_verbose", ["-D"]),
                            ("write_debug", ["--write-dbg"]),
                            ("all_debug", ["-D", "--write-dbg"])]:
            folder, result = run(name, ["--stdout", *extra])
            require(parse(result.stdout) == expected_records, name + " changed or contaminated stdout")
            no_stream_artifacts(folder)

        log_path = root / "requested.log"
        folder, result = run("log_file", ["--stdout", "-v", "--write-dbg", "--log-file", str(log_path)])
        require(parse(result.stdout) == expected_records, "--log-file changed stdout")
        require(log_path.stat().st_size > 0, "explicit log file was not written")
        require(b"[main]" in result.stderr, "--log-file removed stderr diagnostics")
        no_stream_artifacts(folder)

        max_log_path = root / "max-verbose.log"
        folder, result = run("max_log_file", ["--stdout", "-D", "--write-dbg", "--log-file", str(max_log_path)])
        require(parse(result.stdout) == expected_records, "max verbosity/logging contaminated stdout")
        require(max_log_path.stat().st_size > 0, "max-verbose log file was not written")
        no_stream_artifacts(folder)

        fallback_env = dict(env, RAD_NO_PIGZ="1")
        folder, _ = run("zlib_file", local_env=fallback_env)
        require(parse(gzip.decompress((folder / "result.fq.gz").read_bytes())) == expected_records,
                "zlib fallback changed default FASTQ output")

        folder, result = run("empty", ["--stdout", "--min-read-length", "100000"])
        require(result.stdout == b"", "zero accepted reads emitted non-record stdout")
        no_stream_artifacts(folder)

        folder, normal_fasta = run("normal_fasta", source=fasta)
        fasta_records = parse(gzip.decompress((folder / "result.fa.gz").read_bytes()), fasta=True)
        require(len(fasta_records) == len(units), "FASTA fixture did not produce reads")
        folder, result = run("stdout_fasta", ["--stdout", "--write-dbg"], source=fasta)
        require(parse(result.stdout, fasta=True) == fasta_records, "stdout FASTA differs from file output")
        no_stream_artifacts(folder)

        folder, hmm_normal = run("hmm_normal", ["--concat-hmm"])
        require(b"[concat_hmm] Model built:" in hmm_normal.stdout, "HMM fixture did not exercise its writer branch")
        require((folder / "result_boundary_unresolved.fq.gz").exists(), "normal HMM unresolved sidecar changed")
        hmm_records = parse(gzip.decompress((folder / "result.fq.gz").read_bytes()))
        folder, result = run("hmm_stdout", ["--stdout", "--concat-hmm", "--write-dbg"])
        require(b"[concat_hmm] Model built:" in result.stderr, "stdout HMM model did not build")
        require(parse(result.stdout) == hmm_records, "HMM stdout differs from HMM file output")
        no_stream_artifacts(folder)

        folder, result = run("split_conflict", ["--stdout", "--bc-split"], expected=1)
        require(result.stdout == b"" and b"--bc-split" in result.stderr, "stdout/split conflict was not reported cleanly")
        no_stream_artifacts(folder)

        large = root / "many.fastq"
        write_reads(large, (units[i % len(units)] for i in range(2048)))
        folder = root / "early_close"
        args = command(folder, large, ["--stdout", "--chunk-size", "1"])
        with (root / "early_close.stderr").open("wb") as errors:
            process = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=errors, env=env)
            try:
                readable, _, _ = select.select([process.stdout], [], [], 15)
                require(readable, "stream did not flush its first FASTQ record promptly")
                first = b"".join(process.stdout.readline() for _ in range(4))
                require(len(parse(first)) == 1, "early-close fixture did not emit a first FASTQ record")
                process.stdout.close()
                status = process.wait(timeout=15)
            except BaseException:
                process.kill()
                process.wait()
                raise
        require(status == 128 + signal.SIGPIPE,
                "early consumer close did not exit with conventional broken-pipe status: " + str(status))
        no_stream_artifacts(folder)

        if Path("/dev/full").exists():
            folder = root / "write_failure"
            args = command(folder, fastq, ["--stdout"])
            with open("/dev/full", "wb", buffering=0) as full:
                result = subprocess.run(args, stdout=full, stderr=subprocess.PIPE, env=env, timeout=15)
            require(result.returncode == 1, "output failure was ignored, aborted or hung: " + str(result.returncode))
            require(b"[ERROR]" in result.stderr, "output failure was not reported on stderr")
            no_stream_artifacts(folder)

        print("stdout CLI passed: file equivalence, clean diagnostics, debug suppression, HMM, FASTA, tails, empty output and pipe/write failures")


if __name__ == "__main__":
    main()
