# Provenance

GenomeLens was independently conceived and developed by me. No third party has
contributed source code, and the project was not created as coursework,
university employment, sponsored research, or under a grant.

## Data

No human genome, personal genotype file, or personal health record is committed
to this repository, and none has ever been.

Every sample, fixture and figure here is generated from synthetic input, and the
committed economics sample is reproducible end to end:

```bash
python scripts/make_econ_sample.py
```

That builds a synthetic whole-genome VCF, runs the full pipeline against it, and
writes the sample report, the canonical payload and the consistency findings. No
human genome is involved at any step.

## Record

This file is committed and versioned so the statement above travels with the Git
history rather than depending on recollection. The commit log is the primary
creation record.
