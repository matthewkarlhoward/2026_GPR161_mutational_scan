# Releasing and archiving on Zenodo

The repository is connected to Zenodo and **v0.2.0 is deposited**:

| | |
|---|---|
| concept DOI (cite this) | [10.5281/zenodo.23682577](https://doi.org/10.5281/zenodo.23682577) |
| v0.2.0 | [10.5281/zenodo.23754402](https://doi.org/10.5281/zenodo.23754402) |
| v0.1.0 | [10.5281/zenodo.23682578](https://doi.org/10.5281/zenodo.23682578) |

The one-time setup below is done; keep it for the record. Later releases need only the
"Cutting a release" steps.

## One-time setup

1. Sign in to <https://zenodo.org> with the GitHub account that owns this repository.
2. Open the profile menu → **GitHub** (or go straight to
   <https://zenodo.org/account/settings/github/>), press **Sync** if the repository is
   not listed yet, find `matthewkarlhoward/2026_GPR161_mutational_scan` and switch it
   **on**. Zenodo adds a webhook; nothing is archived until a release is published
   *after* this point — earlier tags are not picked up retroactively, though Zenodo
   support can import a past release by hand if you ask.
3. The ORCID is in both files and came through on the record. If you want an
   institutional affiliation shown, add `"affiliation"` inside the `creators` entry of
   `.zenodo.json` (and `affiliation:` under `authors` in `CITATION.cff`) — it is the one
   field still blank, and it will appear from the next release onwards.

## Cutting a release

1. Make sure a fresh clone rebuilds: `pip install -r requirements.txt`, then run the
   pipeline and figure scripts listed in the README.
2. Update `CITATION.cff` — `version:` and `date-released:` — so they agree with the tag
   you are about to push. (`.zenodo.json` carries no version field; Zenodo takes the
   version from the GitHub tag. Add a `"version"` key there only if you want the record
   to say something different from the tag.)
3. Tag and publish a release on GitHub (`v1.0.0`, then `v1.1.0`, `v2.0.0`, …).
   Zenodo archives the tarball of that tag and mints a DOI within a few minutes.

## Versioning afterwards

Zenodo gives the first deposit a **concept DOI** (always resolves to the latest version)
and each release its own **version DOI**. Every later GitHub release is automatically a
new version under the same concept DOI — nothing extra to do, and nothing to re-upload.
Cite the concept DOI in the paper so the citation keeps working as figures are added.

Suggested scheme while the manuscript is in progress:

| tag | when |
|---|---|
| `v0.x.0` | pre-submission snapshots as figures are added |
| `v1.0.0` | the version submitted with the manuscript |
| `v1.x.0` | revisions during review |

## Metadata files

- `CITATION.cff` — GitHub reads it and shows **Cite this repository**. Zenodo reads it
  only when there is no `.zenodo.json`; since this repository has one, the CFF is ignored
  for the deposit and matters for GitHub alone.
- `.zenodo.json` — what Zenodo actually archives with, and the only place for
  `upload_type`, `access_right` and `related_identifiers`. Edit this one to change the
  record; keep the CFF in step by hand so the two do not drift.

The Zenodo record is set to `"license": "cc-by-4.0"`, because most of what is deposited
is data. The code stays MIT inside the archive (`LICENSE`), and the description says so.
If you would rather the record itself read MIT, change that one field.

## Before the submission release, check

- [x] the raw FASTQ are deposited and the SRA accession is written into
      `pipeline/README.md` §1 — BioProject PRJNA1551518, processed 2026-10-27
- [ ] the six `[CONFIRM]` items in `docs/METHODS.md` and `docs/DATA_SOURCES.md` §4 are
      resolved (Dumpling commit, Lilace version and call, and the reporter and
      construct citations)
- [ ] `related_identifiers` in `.zenodo.json` gains the paper's DOI with
      `"relation": "isSupplementTo"` once it exists
- [ ] scores are deposited in MaveDB and that accession is cross-referenced here
