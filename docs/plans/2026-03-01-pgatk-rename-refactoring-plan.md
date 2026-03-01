# pgatk Rename & Refactoring Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Rename `py-pgatk`/`pypgatk` to `pgatk` everywhere, add config auto-loading, fix bugs, and remove dead code.

**Architecture:** 4 independent PRs — PR 1 (rename) merges first, PRs 2-4 are independent of each other. Each PR gets its own branch off `dev`.

**Tech Stack:** Python 3.9+, Click CLI, PyYAML, BioPython, pyOpenMS, MkDocs Material

---

## PR 1: Rename `pypgatk` → `pgatk` (branch: `refactor/rename`)

### Task 1.1: Rename the package directory and CLI module

**Files:**
- Rename: `pypgatk/` → `pgatk/` (entire directory tree)
- Rename: `pgatk/pypgatkc.py` → `pgatk/cli.py`

**Step 1: Rename the directory**

```bash
cd /Users/yperez/work/proteogenomics/pgatk
git mv pypgatk pgatk
```

**Step 2: Rename the CLI module**

```bash
git mv pgatk/pypgatkc.py pgatk/cli.py
```

**Step 3: Verify the rename**

```bash
ls pgatk/
# Should show: __init__.py cli.py cgenomes/ commands/ config/ db/ ensembl/ logs/ proteogenomics/ proteomics/ test_all.bed test_annotated.vcf testdata/ tests/ toolbox/
```

---

### Task 1.2: Update `__init__.py`

**Files:**
- Modify: `pgatk/__init__.py`

**Step 1: Update the package name**

Change `name = "pypgatk"` to `name = "pgatk"`.

---

### Task 1.3: Update `cli.py` (formerly `pypgatkc.py`)

**Files:**
- Modify: `pgatk/cli.py`

**Step 1: Update all imports from `pypgatk` to `pgatk`**

Replace all 14 import lines:
```python
from pgatk.commands import ensembl_downloader as ensembl_downloader_cmd
from pgatk.commands import ensembl_database as ensembl_database_cmd
from pgatk.commands import cosmic_downloader as cosmic_downloader_cmd
from pgatk.commands import cbioportal_downloader as cbioportal_downloader_cmd
from pgatk.commands import cosmic_to_proteindb as cosmic_to_proteindb_cmd
from pgatk.commands import cbioportal_to_proteindb as cbioportal_to_proteindb_cmd
from pgatk.commands import threeframe_translation as threeframe_translation_cmd
from pgatk.commands import vcf_to_proteindb as vcf_to_proteindb_cmd
from pgatk.commands import dnaseq_to_proteindb as dnase_to_proteindb_cmd
from pgatk.commands import proteindb_decoy as proteindb_decoy_cmd
from pgatk.commands import peptide_class_fdr as peptide_class_fdr_cmd
from pgatk.commands import validate_peptides as validate_peptides_cmd
from pgatk.commands import mztab_class_fdr as mztab_class_fdr_cmd
from pgatk.commands import blast_get_position as blast_get_position_cmd
```

**Step 2: Update docstring**

Change `"py-pgatk"` references to `"pgatk"` in the module docstring and CLI help text.

---

### Task 1.4: Update all Python imports across all modules

**Files (all files under `pgatk/` that import from `pypgatk`):**

Use a global find-and-replace: `from pypgatk.` → `from pgatk.`

These files need updating:
- `pgatk/commands/blast_get_position.py`
- `pgatk/commands/cbioportal_downloader.py`
- `pgatk/commands/cbioportal_to_proteindb.py`
- `pgatk/commands/cosmic_downloader.py`
- `pgatk/commands/cosmic_to_proteindb.py`
- `pgatk/commands/dnaseq_to_proteindb.py`
- `pgatk/commands/ensembl_database.py`
- `pgatk/commands/ensembl_downloader.py`
- `pgatk/commands/mztab_class_fdr.py`
- `pgatk/commands/peptide_class_fdr.py`
- `pgatk/commands/proteindb_decoy.py`
- `pgatk/commands/threeframe_translation.py`
- `pgatk/commands/validate_peptides.py`
- `pgatk/commands/vcf_to_proteindb.py`
- `pgatk/cgenomes/cbioportal_downloader.py`
- `pgatk/cgenomes/cgenomes_proteindb.py`
- `pgatk/cgenomes/cosmic_downloader.py`
- `pgatk/ensembl/data_downloader.py`
- `pgatk/ensembl/ensembl.py`
- `pgatk/ensembl/exceptions.py`
- `pgatk/proteogenomics/blast_get_position.py`
- `pgatk/proteogenomics/mztab_class_fdr.py`
- `pgatk/proteogenomics/spectrumai.py`
- `pgatk/proteomics/openms.py`
- `pgatk/proteomics/db/protein_database_decoy.py`
- `pgatk/toolbox/general.py`

**Step 1: Run a sed-like replacement across all `.py` files**

In every `.py` file under `pgatk/`, replace `from pypgatk.` with `from pgatk.` and `import pypgatk` with `import pgatk`.

**Step 2: Verify no remaining references**

```bash
grep -r "pypgatk" pgatk/ --include="*.py"
# Should return nothing (except possibly string literals in help text)
```

---

### Task 1.5: Update test files

**Files:**
- Modify: `pgatk/tests/pypgatk_tests.py`
- Modify: `pgatk/tests/classes_tests.py`

**Step 1: Update test imports**

In `pypgatk_tests.py`:
```python
# Before
from pypgatk.pypgatkc import cli
# After
from pgatk.cli import cli
```

In `classes_tests.py`:
```python
# Before
from pypgatk.proteomics.models import PYGPATK_ENZYMES
# After
from pgatk.proteomics.models import PYGPATK_ENZYMES
```

**Step 2: Rename the constant `PYGPATK_ENZYMES` → `PGATK_ENZYMES`**

In `pgatk/proteomics/models.py`, rename the constant:
```python
# Before
PYGPATK_ENZYMES = ...
# After
PGATK_ENZYMES = ...
```

Update all references (in `models.py`, `classes_tests.py`, `protein_database_decoy.py`, and any command files).

---

### Task 1.6: Update `pyproject.toml`

**Files:**
- Modify: `pyproject.toml`

**Step 1: Update all references**

```toml
[project]
name = "pgatk"

[project.scripts]
pgatk = "pgatk.cli:main"

[project.urls]
Homepage = "http://github.com/bigbio/pgatk"
Repository = "http://github.com/bigbio/pgatk"
Issues = "http://github.com/bigbio/pgatk/issues"

[tool.setuptools.packages.find]
include = ["pgatk*"]

[tool.setuptools.package-data]
pgatk = ["config/*.yaml", "config/*.json"]
```

---

### Task 1.7: Update non-Python files

**Files:**
- Modify: `MANIFEST.in`
- Modify: `Dockerfile`
- Modify: `.gitignore`
- Modify: `conda-enviroment.yaml`

**Step 1: Update MANIFEST.in**

```
recursive-include pgatk *.json
recursive-include pgatk *.yaml
```

**Step 2: Update Dockerfile**

Replace `pypgatk` references with `pgatk`:
- `LABEL software="pgatk"`, `container="pgatk"`
- `ENV PATH=$PATH:/tool/source/pgatk/`
- `RUN chmod +x /tool/source/pgatk/cli.py`
- Update git clone URL to `bigbio/pgatk`

**Step 3: Update .gitignore**

Replace `pypgatk/testdata/` paths with `pgatk/testdata/`.

**Step 4: Update conda-enviroment.yaml**

Change `name: py-pgatk` to `name: pgatk`.

---

### Task 1.8: Update GitHub Actions workflows

**Files:**
- Modify: `.github/workflows/pythonapp.yml`
- Modify: `.github/workflows/pythonpackage.yml`

**Step 1: Update test directory references**

In both files, change:
```yaml
# Before
cd pypgatk
# After
cd pgatk
```

---

### Task 1.9: Update documentation files

**Files:**
- Modify: `mkdocs.yml`
- Modify: `README.md`
- Modify: `docs/index.md`
- Modify: `docs/installation.md`
- Modify: `docs/pypgatk.md` (rename to `docs/pgatk-cli.md`)
- Modify: `docs/changelog.md`

**Step 1: Update mkdocs.yml**

```yaml
repo_name: bigbio/pgatk
repo_url: https://github.com/bigbio/pgatk
# Update nav
nav:
  - Home: index.md
  - Introduction: introduction.md
  - Installation: installation.md
  - pgatk CLI: pgatk-cli.md
  - File Formats: formats.md
  - Changelog: changelog.md
  - Support: support.md
# Update social links
- icon: fontawesome/brands/github
  link: https://github.com/bigbio/pgatk
- icon: fontawesome/brands/python
  link: https://pypi.org/project/pgatk/
```

**Step 2: Rename docs/pypgatk.md → docs/pgatk-cli.md**

```bash
git mv docs/pypgatk.md docs/pgatk-cli.md
```

**Step 3: Update all doc content**

In all doc files, replace:
- `pypgatk` (the CLI command) → `pgatk`
- `py-pgatk` (repo name) → `pgatk`
- `pip install pypgatk` → `pip install pgatk`
- Badge URLs: `bigbio/py-pgatk` → `bigbio/pgatk`
- PyPI badges: `pypgatk` → `pgatk`

---

### Task 1.10: Run tests and verify

**Step 1: Install in editable mode**

```bash
cd /Users/yperez/work/proteogenomics/pgatk
pip install -e .
```

**Step 2: Run tests**

```bash
cd pgatk
python -m unittest discover -s tests -p "*_tests.py" -v
```

Expected: All tests pass.

**Step 3: Verify CLI entry point**

```bash
pgatk --version
# Expected: pgatk, version 0.0.26
```

**Step 4: Commit**

```bash
git add -A
git commit -m "refactor: rename pypgatk to pgatk

- Rename package directory pypgatk/ → pgatk/
- Rename CLI module pypgatkc.py → cli.py
- Update all Python imports from pypgatk to pgatk
- Rename constant PYGPATK_ENZYMES → PGATK_ENZYMES
- Update pyproject.toml (name, entry point, URLs, package config)
- Update Dockerfile, MANIFEST.in, .gitignore, conda env
- Update GitHub Actions workflows
- Update MkDocs documentation
- CLI command: pypgatk → pgatk
- PyPI package: pypgatk → pgatk"
```

---

## PR 2: Config Registry (branch: `refactor/config-registry`)

### Task 2.1: Create `pgatk/config/registry.py`

**Files:**
- Create: `pgatk/config/registry.py`

**Step 1: Write the registry module**

```python
"""
Central config registry: maps config keys to bundled default YAML files.
Enables auto-loading of defaults when user doesn't pass -c.
"""
import os

from pgatk.toolbox.general import read_yaml_from_file

CONFIG_DIR = os.path.dirname(__file__)

COMMAND_CONFIGS = {
    "ensembl_downloader": "ensembl_downloader_config.yaml",
    "ensembl_config": "ensembl_config.yaml",
    "cosmic": "cosmic_config.yaml",
    "cbioportal": "cbioportal_config.yaml",
    "protein_decoy": "protein_decoy.yaml",
    "openms_analysis": "openms_analysis.yaml",
}


def load_config(config_key, user_config_file=None):
    """
    Load configuration for a command.

    If user_config_file is provided, loads that file.
    Otherwise, loads the bundled default config for the given key.

    :param config_key: key in COMMAND_CONFIGS
    :param user_config_file: optional path to user-provided config file
    :return: parsed YAML config dict, or empty dict if no config found
    """
    if user_config_file is not None:
        return read_yaml_from_file(user_config_file)

    default_file = COMMAND_CONFIGS.get(config_key)
    if default_file is None:
        return {}
    return read_yaml_from_file(os.path.join(CONFIG_DIR, default_file))
```

---

### Task 2.2: Update command files to use `load_config`

**Files (all 14 command files in `pgatk/commands/`):**

For each command file, replace:

```python
# Before
from pgatk.toolbox.general import read_yaml_from_file

# ... in the command function:
config_data = None
if config_file is not None:
    config_data = read_yaml_from_file(config_file)
```

With:

```python
# After
from pgatk.config.registry import load_config

# ... in the command function:
config_data = load_config("<config_key>", config_file)
```

**Config key mapping per command file:**

| Command file | Config key |
|---|---|
| `ensembl_downloader.py` | `"ensembl_downloader"` |
| `ensembl_database.py` | `"ensembl_config"` |
| `vcf_to_proteindb.py` | `"ensembl_config"` |
| `dnaseq_to_proteindb.py` | `"ensembl_config"` |
| `threeframe_translation.py` | `"ensembl_config"` |
| `cosmic_downloader.py` | `"cosmic"` |
| `cosmic_to_proteindb.py` | `"cosmic"` |
| `cbioportal_downloader.py` | `"cbioportal"` |
| `cbioportal_to_proteindb.py` | `"cbioportal"` |
| `proteindb_decoy.py` | `"protein_decoy"` |
| `peptide_class_fdr.py` | `"openms_analysis"` |
| `validate_peptides.py` | `"ensembl_config"` |
| `mztab_class_fdr.py` | `"ensembl_config"` |
| `blast_get_position.py` | (no config needed — uses None) |

---

### Task 2.3: Update tests that pass `--config_file`

**Files:**
- Modify: `pgatk/tests/pypgatk_tests.py`

Some tests currently pass `--config_file` explicitly. These should still work (user override).
Add a new test that runs a command without `--config_file` to verify auto-loading:

```python
def test_generate_decoy_database_autoconfig(self):
    """Test that generate-decoy works without explicit --config_file (auto-loads bundled default)."""
    runner = CliRunner()
    result = runner.invoke(cli,
                           ['generate-decoy', '-in', 'testdata/test_db.fa',
                            '-out', 'testdata/output_decoy.fa',
                            '--method', 'protein-reverse'])
    self.assertEqual(result.exit_code, 0)
```

Note: `test_generate_decoy_database_noconfig` already tests this case — verify it still passes.

**Step 1: Run tests**

```bash
cd pgatk
python -m unittest discover -s tests -p "*_tests.py" -v
```

**Step 2: Commit**

```bash
git add pgatk/config/registry.py pgatk/commands/
git commit -m "refactor: add config registry for auto-loading bundled defaults

- Add pgatk/config/registry.py with central COMMAND_CONFIGS mapping
- Update all 14 command files to use load_config()
- -c flag is now optional: bundled defaults are loaded automatically
- User can still override with -c path/to/custom_config.yaml"
```

---

## PR 3: Bug Fixes (branch: `fix/bugs`)

### Task 3.1: Fix `_filter_write_idxml_with_df` writing unfiltered data

**Files:**
- Modify: `pgatk/proteomics/openms.py:437`

**Step 1: Fix the store call**

```python
# Before (line 437)
idxml_parser().store(output_file, prot_ids, pep_ids)
# After
idxml_parser().store(output_file, prot_ids, new_pep_ids)
```

---

### Task 3.2: Fix ISO-8859 encoding typo

**Files:**
- Modify: `pgatk/cgenomes/cosmic_downloader.py:156`

**Step 1: Fix the encoding name**

The en-dash (U+2013) in `'ISO-8859–1'` must be replaced with a regular hyphen:

```python
# Before (line 156)
outfile.write(gzip.decompress(open(local_file, 'rb').read()).decode('ISO-8859–1'))
# After
outfile.write(gzip.decompress(open(local_file, 'rb').read()).decode('ISO-8859-1'))
```

---

### Task 3.3: Fix `download_file` log inversion bug

**Files:**
- Modify: `pgatk/toolbox/general.py:192-193`

**Step 1: Fix the condition**

```python
# Before (line 192-193)
if log is not None:
    log = logging
# After
if log is None:
    log = logging
```

---

### Task 3.4: Fix enzyme name typos in models.py

**Files:**
- Modify: `pgatk/proteomics/models.py:21-22`

**Step 1: Fix the typos**

```python
# Before (lines 21-22 — "cleavege" is misspelled)
'Trypsin/P': {'cleavage rule': '(?<=[KRX])', 'PSIID': 'MS:1001313', 'cleavege sites': 'KR'},
'V8-DE': {'cleavage rule': '(?<=[DBEZX])(?!P)', 'PSIID': 'MS:1001314', 'cleavege sites': 'DBEZX'},
# After
'Trypsin/P': {'cleavage rule': '(?<=[KRX])', 'PSIID': 'MS:1001313', 'cleavage sites': 'KR'},
'V8-DE': {'cleavage rule': '(?<=[DBEZX])(?!P)', 'PSIID': 'MS:1001314', 'cleavage sites': 'DBEZX'},
```

**Important:** Check if `'cleavege sites'` is used as a dictionary key anywhere — if so, update those references too. Search for `cleavege` across the codebase.

---

### Task 3.5: Replace `print_help()` with Click's `required=True`

**Files (11 command files):**

For each command file listed below, add `required=True` to the relevant `@click.option()` decorators and remove the `print_help()` check.

| File | Options to make required |
|---|---|
| `ensembl_database.py` | `--input_fasta` |
| `vcf_to_proteindb.py` | `--input_fasta`, `--vcf`, `--gene_annotations_gtf` |
| `cbioportal_to_proteindb.py` | `--input_mutation`, `--input_cds`, `--output_db` |
| `cosmic_to_proteindb.py` | `--input_mutation`, `--input_genes`, `--output_db` |
| `threeframe_translation.py` | `--input_fasta` |
| `proteindb_decoy.py` | `--input_database` |
| `peptide_class_fdr.py` | `--input_file`, `--output_file` |
| `dnaseq_to_proteindb.py` | `--input_fasta` |
| `blast_get_position.py` | `--input_psm_to_blast`, `--input_reference_database`, `--output_psm` |
| `mztab_class_fdr.py` | `--input_mztab`, `--outfile_name` |

**Special case: `ensembl_downloader.py`** — This command requires `--taxonomy` OR `--ensembl_name` (either one). Can't use `required=True` directly. Use a Click callback or keep a validation check (but raise `click.UsageError` instead of `print_help()`):

```python
if taxonomy is None and ensembl_name is None:
    raise click.UsageError("Either --taxonomy or --ensembl_name is required.")
```

**Special case: `validate_peptides.py`** — Review the specific validation logic before converting.

**Step 1: For each file, add `required=True` and remove `print_help()` import and call**

Example pattern:
```python
# Before
from pgatk.commands.utils import print_help

@click.option('-in', '--input_fasta', help='...')
def command(ctx, config_file, input_fasta, ...):
    if input_fasta is None:
        print_help()

# After (remove print_help import, add required=True)
@click.option('-in', '--input_fasta', help='...', required=True)
def command(ctx, config_file, input_fasta, ...):
    # print_help() call removed
```

**Step 2: Remove `print_help` from `pgatk/commands/utils.py`** if no longer used anywhere.

**Step 3: Run tests**

```bash
cd pgatk
python -m unittest discover -s tests -p "*_tests.py" -v
```

---

### Task 3.6: Replace deprecated `Bio.pairwise2`

**Files:**
- Modify: `pgatk/proteogenomics/blast_get_position.py`

**Step 1: Replace import**

```python
# Before
from Bio import pairwise2, SeqIO
# After
from Bio import SeqIO
from Bio.Align import PairwiseAligner
```

**Step 2: Replace usage in `peptide_blast_protein` function**

```python
# Before
score = pairwise2.align.localms(sequenceA=fasta, sequenceB=peptide,
                                match=1, mismatch=0, open=-2, extend=-2, score_only=True)
if score == length-1:
    alignment = pairwise2.align.localms(sequenceA=fasta, sequenceB=peptide,
                                        match=1, mismatch=0, open=-2, extend=-2)[0]

# After
aligner = PairwiseAligner()
aligner.mode = 'local'
aligner.match_score = 1
aligner.mismatch_score = 0
aligner.open_gap_score = -2
aligner.extend_gap_score = -2
alignments = aligner.align(fasta, peptide)
score = alignments.score if alignments else 0
if score == length - 1:
    alignment = alignments[0]
    # Update attribute access: alignment.aligned gives coordinate pairs
    # Need to adapt get_details() calls to work with new alignment format
```

**Note:** The `PairwiseAligner` API differs from `pairwise2`. The alignment object attributes change:
- `alignment.seqA` / `alignment.seqB` → access via `alignment.target` / `alignment.query`
- `alignment.start` / `alignment.end` → access via `alignment.aligned` coordinate arrays

Carefully test the `get_details()` function calls to ensure they still work with the new alignment format.

**Step 3: Commit**

```bash
git add -A
git commit -m "fix: resolve multiple bugs

- Fix _filter_write_idxml_with_df writing unfiltered data (openms.py)
- Fix ISO-8859 encoding en-dash typo (cosmic_downloader.py)
- Fix download_file log inversion bug (general.py)
- Fix enzyme name typos 'cleavege' → 'cleavage' (models.py)
- Replace print_help() with Click required=True (11 command files)
- Replace deprecated Bio.pairwise2 with PairwiseAligner (blast_get_position.py)"
```

---

## PR 4: Cleanup (branch: `refactor/cleanup`)

### Task 4.1: Remove dead `Species` and `SpeciesService` classes

**Files:**
- Modify: `pgatk/ensembl/models.py`

**Step 1: Verify they are unused**

```bash
grep -r "Species\|SpeciesService" pgatk/ --include="*.py" -l
# Should only show pgatk/ensembl/models.py
```

**Step 2: Remove both classes** (lines 9-180 approximately). Keep any other code in the file (e.g., imports used elsewhere).

---

### Task 4.2: Remove dead code after return statement

**Files:**
- Modify: `pgatk/toolbox/general.py`

**Step 1: Delete commented-out code** at lines 241-253 (after `return downloaded_file` in `download_file()`).

---

### Task 4.3: Dynamic version from package metadata

**Files:**
- Modify: `pgatk/cli.py`

**Step 1: Use importlib.metadata**

```python
# Before
@click.version_option(version='0.0.26')
def cli():

# After
from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("pgatk")
except PackageNotFoundError:
    __version__ = "dev"

@click.version_option(version=__version__)
def cli():
```

---

### Task 4.4: Fix CWD log file creation

**Files:**
- Modify: `pgatk/toolbox/general.py`

**Step 1: Replace FileHandler with StreamHandler**

In `ParameterConfiguration.__init__`, replace the file-based logging with stderr-based logging:

```python
# Before (lines 72-89)
self._log_handlers = []
log_handlers_prefix = self._ROOT_CONFIG_NAME + '-'
log_handlers_extension = '.log'

self._logger = logging.getLogger(__name__)
self._logger.setLevel(getattr(logging, self._log_level))
self._log_files = []
for llevel, lformat in self._logger_formatters.items():
    logfile = os.path.join(log_handlers_prefix + llevel.lower() + log_handlers_extension)
    lformatter = logging.Formatter(lformat)
    lhandler = logging.FileHandler(logfile, mode='w')
    lhandler.setLevel(getattr(logging, llevel))
    lhandler.setFormatter(lformatter)
    self._log_handlers.append(lhandler)
    self._logger.addHandler(lhandler)
    self._log_files.append(logfile)

# After
self._log_handlers = []
self._logger = logging.getLogger(__name__)
self._logger.setLevel(getattr(logging, self._log_level))
self._log_files = []

# Log to stderr instead of creating files in CWD
for llevel, lformat in self._logger_formatters.items():
    lformatter = logging.Formatter(lformat)
    lhandler = logging.StreamHandler()
    lhandler.setLevel(getattr(logging, llevel))
    lhandler.setFormatter(lformatter)
    self._log_handlers.append(lhandler)
    self._logger.addHandler(lhandler)
```

**Step 2: Clean up `get_session_log_files()`**

This method returns `self._log_files` which will now be empty. Keep the method (in case downstream code calls it) but it will return `[]`.

---

### Task 4.5: Run full test suite and commit

**Step 1: Run tests**

```bash
cd pgatk
python -m unittest discover -s tests -p "*_tests.py" -v
```

**Step 2: Commit**

```bash
git add -A
git commit -m "refactor: remove dead code and improve logging

- Remove unused Species and SpeciesService classes (ensembl/models.py)
- Remove commented-out dead code after return (toolbox/general.py)
- Use dynamic version from package metadata instead of hardcoded '0.0.26'
- Log to stderr instead of creating log files in CWD"
```

---

## Summary: Branch and PR Strategy

| PR | Branch | Base | Description |
|---|---|---|---|
| 1 | `refactor/rename` | `dev` | Full rename: pypgatk → pgatk |
| 2 | `refactor/config-registry` | `refactor/rename` | Auto-load bundled config defaults |
| 3 | `fix/bugs` | `refactor/rename` | Fix 6 bugs |
| 4 | `refactor/cleanup` | `refactor/rename` | Remove dead code, fix logging |

PR 1 merges first → then PRs 2, 3, 4 can merge in any order.

## Post-Merge (Manual)

1. Rename GitHub repo: Settings → Repository name → `pgatk`
2. Publish new `pgatk` package to PyPI
3. Update Bioconda recipe
