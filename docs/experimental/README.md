# Experimental investigations

These notes preserve completed or ongoing engineering investigations. They are
evidence records, not supported RAD features: code may remain on a separate
branch or in a private dependency fork, and a positive microbenchmark does not
mean the approach has been adopted.

| Experiment | Dates | Status | Outcome |
| --- | --- | --- | --- |
| [ISA-L igzip 2.31.1 on Apple Silicon](igzip-2.31.1-macos.md) | 2026-07-26–27 | Concluded; retained privately | Patched compression was 45.9% faster in isolation, but the best RAD configuration improved median end-to-end time by only 1.9%; pigz remains the default. |
