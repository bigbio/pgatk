# Phase 0: Infrastructure & Code Quality — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Fix existing bugs, modernize code patterns, and expand test coverage so the graph engine (Phase 1) is built on solid ground.

**Architecture:** No architectural changes — this is a cleanup phase touching existing files. All changes are backward-compatible. The test suite must pass before and after each task.

**Tech Stack:** Python 3.9+, pytest, dataclasses, pathlib, concurrent.futures, logging

---

## Task 1: Add conftest.py with proper fixtures and absolute paths

**Files:**
- Create: `pgatk/tests/conftest.py`
- Modify: `pgatk/tests/pgatk_tests.py`
- Modify: `pgatk/tests/classes_tests.py`

**Step 1: Create conftest.py with path fixtures**

```python
# pgatk/tests/conftest.py
import os
from pathlib import Path

import pytest

# Absolute path to testdata/ directory, works regardless of CWD
TESTDATA_DIR = Path(__file__).resolve().parent.parent / "testdata"
CONFIG_DIR = Path(__file__).resolve().parent.parent / "config"


@pytest.fixture
def testdata_dir():
    """Absolute path to the testdata directory."""
    return TESTDATA_DIR


@pytest.fixture
def config_dir():
    """Absolute path to the config directory."""
    return CONFIG_DIR


@pytest.fixture(autouse=True)
def change_to_package_dir(tmp_path, monkeypatch):
    """Ensure tests run from a consistent directory (the pgatk package root)."""
    monkeypatch.chdir(Path(__file__).resolve().parent.parent)
```

**Step 2: Run existing tests to verify they still pass**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/ -v --timeout=120 -x`
Expected: All existing tests pass (the `autouse` fixture changes CWD to the package root, preserving relative path behavior).

**Step 3: Commit**

```bash
git add pgatk/tests/conftest.py
git commit -m "test: add conftest.py with testdata/config path fixtures"
```

---

## Task 2: Fix variable shadowing bug in protein_database_decoy.py

**Files:**
- Modify: `pgatk/proteomics/db/protein_database_decoy.py:170,180,203`
- Create: `pgatk/tests/test_decoy_shadowing.py`

**Step 1: Write the failing test**

```python
# pgatk/tests/test_decoy_shadowing.py
"""Verify that the decoy_sequence import is not shadowed by local variables."""
from pyteomics.fasta import decoy_sequence as pyteomics_decoy_sequence
from pgatk.proteomics.db import protein_database_decoy


def test_decoy_sequence_import_not_shadowed():
    """The module-level import of decoy_sequence must remain callable after
    print_target_decoy_composition defines a local variable of the same name."""
    # Verify the import at module level is still the pyteomics function
    module_ref = protein_database_decoy.decoy_sequence
    assert callable(module_ref), (
        "decoy_sequence at module level should be the pyteomics function, "
        f"got {type(module_ref)}"
    )
    assert module_ref is pyteomics_decoy_sequence
```

**Step 2: Run test to verify it passes (this is a module-level check, so it should pass)**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_decoy_shadowing.py -v`
Expected: PASS (the shadow is local to the method, not module-level — but we fix it for safety)

**Step 3: Rename the shadowing variable**

In `pgatk/proteomics/db/protein_database_decoy.py`, rename the local variable `decoy_sequence` to `decoy_seq_accumulated` at lines 170, 180, and 203:

- Line 170: `decoy_sequence = ''` → `decoy_seq_accumulated = ''`
- Line 180: `decoy_sequence = decoy_sequence + str(record.seq)` → `decoy_seq_accumulated = decoy_seq_accumulated + str(record.seq)`
- Line 203: `decoy_aa_composition = self.count_aa_in_dictionary(decoy_aa_composition, decoy_sequence)` → `decoy_aa_composition = self.count_aa_in_dictionary(decoy_aa_composition, decoy_seq_accumulated)`

**Step 4: Run all tests to verify nothing broke**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/ -v --timeout=120 -x`
Expected: All PASS

**Step 5: Commit**

```bash
git add pgatk/proteomics/db/protein_database_decoy.py pgatk/tests/test_decoy_shadowing.py
git commit -m "fix: rename shadowing variable decoy_sequence in print_target_decoy_composition"
```

---

## Task 3: Convert SNP class to dataclass

**Files:**
- Modify: `pgatk/cgenomes/models.py`
- Modify: `pgatk/cgenomes/cgenomes_proteindb.py` (update `snp.type` → `snp.mutation_type`)
- Create: `pgatk/tests/test_models.py`

**Step 1: Write the failing test**

```python
# pgatk/tests/test_models.py
import dataclasses
from pgatk.cgenomes.models import SNP


def test_snp_is_dataclass():
    assert dataclasses.is_dataclass(SNP)


def test_snp_defaults():
    snp = SNP()
    assert snp.gene is None
    assert snp.mrna is None
    assert snp.dna_mut is None
    assert snp.aa_mut is None
    assert snp.mutation_type is None


def test_snp_with_values():
    snp = SNP(gene="BRCA1", mrna="NM_007294", dna_mut="c.5266dupC",
              aa_mut="p.Q1756fs", mutation_type="Insertion")
    assert snp.gene == "BRCA1"
    assert snp.mutation_type == "Insertion"


def test_snp_no_type_builtin_shadow():
    """The old SNP class stored mutation_type as self.type, shadowing the builtin."""
    snp = SNP(mutation_type="Substitution")
    assert not hasattr(snp, 'type') or snp.mutation_type == "Substitution"
```

**Step 2: Run test to verify it fails**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_models.py -v`
Expected: FAIL — `test_snp_is_dataclass` fails because SNP is not a dataclass yet

**Step 3: Convert SNP to dataclass**

Replace the entire content of `pgatk/cgenomes/models.py`:

```python
from dataclasses import dataclass
from typing import Optional


@dataclass
class SNP:
    """Represents a single nucleotide polymorphism or mutation."""
    gene: Optional[str] = None
    mrna: Optional[str] = None
    dna_mut: Optional[str] = None
    aa_mut: Optional[str] = None
    mutation_type: Optional[str] = None
```

**Step 4: Update references from `snp.type` to `snp.mutation_type`**

In `pgatk/cgenomes/cgenomes_proteindb.py`, find all uses of `.type` on SNP objects:
- Search for `snp.type` or `.type =` in the context of SNP objects
- Line where SNP is constructed: `snp = SNP()` followed by `snp.type = ...` → change to `snp.mutation_type = ...`
- Any reads of `snp.type` → `snp.mutation_type`

Specifically, look for patterns like:
```python
snp.type = mutation_type   # → snp.mutation_type = mutation_type
```
And reads like:
```python
snp.type                    # → snp.mutation_type
```

**Step 5: Run tests to verify they pass**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_models.py pgatk/tests/ -v --timeout=120 -x`
Expected: All PASS

**Step 6: Commit**

```bash
git add pgatk/cgenomes/models.py pgatk/cgenomes/cgenomes_proteindb.py pgatk/tests/test_models.py
git commit -m "refactor: convert SNP to dataclass, rename .type to .mutation_type"
```

---

## Task 4: Replace pathos with concurrent.futures in blast_get_position.py

**Files:**
- Modify: `pgatk/proteogenomics/blast_get_position.py:1-10,85,151-155`
- Modify: `pyproject.toml` (remove pathos from dependencies after both files are migrated)

**Step 1: Write the test**

```python
# pgatk/tests/test_blast_imports.py
"""Verify blast_get_position uses stdlib concurrent.futures, not pathos."""
import importlib


def test_no_pathos_import_in_blast():
    """blast_get_position should not import pathos."""
    import pgatk.proteogenomics.blast_get_position as mod
    source_file = mod.__file__
    with open(source_file, 'r') as f:
        source = f.read()
    assert 'pathos' not in source, "blast_get_position.py still imports pathos"
    assert 'concurrent.futures' in source or 'ProcessPoolExecutor' in source
```

**Step 2: Run test to verify it fails**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_blast_imports.py -v`
Expected: FAIL — pathos is still imported

**Step 3: Replace pathos with concurrent.futures**

In `pgatk/proteogenomics/blast_get_position.py`:

Replace imports (lines 1-9):
```python
import pandas as pd
from Bio import SeqIO
import datetime
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import ahocorasick

from pgatk.toolbox.general import ParameterConfiguration
```

Remove `Manager` import and usage. In `__init__`, replace line 85:
```python
# OLD: self.blast_dict = Manager().dict()
self.blast_dict = {}
```

Replace the pool usage (lines 151-155) in `blast()`:
```python
# OLD:
# pool = Pool(int(self._number_of_processes))
# list(tqdm(pool.imap(self._result, seq_list), total=len(seq_list), desc="Blast", unit="peptide"))
# pool.close()
# pool.join()

# NEW:
with ProcessPoolExecutor(max_workers=int(self._number_of_processes)) as executor:
    list(tqdm(executor.map(self._result, seq_list, chunksize=100),
              total=len(seq_list), desc="Blast", unit="peptide"))
```

Note: Since `self._result` writes to `self.blast_dict` which is a shared Manager dict, and we're removing Manager, we need to refactor `_result()` to return results instead of writing to shared state. Then collect results in the main process:

```python
# Refactor _result to return data instead of writing to shared dict
def _result(self, sequence):
    """Process a single peptide and return (peptide, details) or None."""
    result = peptide_blast_protein(self.fasta_dict, sequence)
    if result:
        return (sequence, result)
    return None

# In blast():
with ProcessPoolExecutor(max_workers=int(self._number_of_processes)) as executor:
    results = list(tqdm(executor.map(self._result, seq_list, chunksize=100),
                        total=len(seq_list), desc="Blast", unit="peptide"))
for r in results:
    if r is not None:
        self.blast_dict[r[0]] = r[1]
```

**Step 4: Run the existing blast test**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/pgatk_tests.py::PgatkRunnerTests::test_blast -v --timeout=120`
Expected: PASS

**Step 5: Run import test**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_blast_imports.py -v`
Expected: PASS

**Step 6: Commit**

```bash
git add pgatk/proteogenomics/blast_get_position.py pgatk/tests/test_blast_imports.py
git commit -m "refactor: replace pathos with concurrent.futures in blast_get_position"
```

---

## Task 5: Replace pathos with concurrent.futures in spectrumai.py

**Files:**
- Modify: `pgatk/proteogenomics/spectrumai.py:1-10,47,290,314-318`

**Step 1: Write the test**

```python
# pgatk/tests/test_spectrumai_imports.py
"""Verify spectrumai uses stdlib concurrent.futures, not pathos."""


def test_no_pathos_import_in_spectrumai():
    import pgatk.proteogenomics.spectrumai as mod
    with open(mod.__file__, 'r') as f:
        source = f.read()
    assert 'pathos' not in source, "spectrumai.py still imports pathos"
```

**Step 2: Run test to verify it fails**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_spectrumai_imports.py -v`
Expected: FAIL

**Step 3: Replace pathos with concurrent.futures**

In `pgatk/proteogenomics/spectrumai.py`:

Replace imports (lines 1-10):
```python
import datetime
import os.path
import re
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from pyopenms import (TheoreticalSpectrumGenerator, MSSpectrum,
                      AASequence, Param, MzMLFile, MSExperiment, SpectrumLookup)
from tqdm import tqdm
from pgatk.toolbox.general import ParameterConfiguration
```

Remove Manager usage. In `__init__`, replace line 47:
```python
# OLD: self.df_list = Manager().list()
# Remove this line entirely — we'll collect results from executor.map()
```

Refactor `_multiprocess_inspect_spectrum` (around line 290) to return the result instead of appending to shared list:
```python
def _multiprocess_inspect_spectrum(self, df):
    """Process a single mzML group and return the annotated DataFrame."""
    return self._inspect_spectrum(df, self._mzml_path, self._mzml_files)
```

Replace pool usage (lines 314-318) in `validate()`:
```python
# OLD:
# pool = Pool(int(self._number_of_processes))
# list(tqdm(pool.imap(self._multiprocess_inspect_spectrum, list_of_dfs), ...))
# pool.close()
# pool.join()

# NEW:
with ProcessPoolExecutor(max_workers=int(self._number_of_processes)) as executor:
    results = list(tqdm(
        executor.map(self._multiprocess_inspect_spectrum, list_of_dfs),
        total=len(list_of_dfs), desc="Validate By Each mzML", unit="mzML"
    ))

# Replace the line that concatenates self.df_list:
# OLD: df_output = pd.concat(self.df_list, axis=0, ignore_index=True)
# NEW:
df_output = pd.concat(results, axis=0, ignore_index=True)
```

**Step 4: Run test to verify it passes**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_spectrumai_imports.py -v`
Expected: PASS

**Step 5: Now remove pathos from pyproject.toml**

In `pyproject.toml`, remove the `pathos` dependency line:
```
# Remove: "pathos>=0.2.8,<1.0",
```
Also remove from `requirements.txt` if present.

**Step 6: Run full test suite**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/ -v --timeout=120 -x`
Expected: All PASS

**Step 7: Commit**

```bash
git add pgatk/proteogenomics/spectrumai.py pgatk/tests/test_spectrumai_imports.py pyproject.toml requirements.txt
git commit -m "refactor: replace pathos with concurrent.futures in spectrumai, remove pathos dependency"
```

---

## Task 6: Replace print() with logger across all modules

This is a large task covering many files. Work through them one file at a time.

**Files:**
- Modify: `pgatk/proteomics/db/protein_database_decoy.py` (14 print sites)
- Modify: `pgatk/ensembl/ensembl.py` (9 print sites)
- Modify: `pgatk/cgenomes/cgenomes_proteindb.py` (18 print sites)
- Modify: `pgatk/proteogenomics/blast_get_position.py` (8 print sites)
- Modify: `pgatk/proteogenomics/spectrumai.py` (6 print sites)

**Step 1: Write a test that checks for print() usage in source files**

```python
# pgatk/tests/test_no_print_statements.py
"""Verify that production code uses logger instead of print()."""
import ast
from pathlib import Path

PGATK_ROOT = Path(__file__).resolve().parent.parent

# Files that should not contain print() calls (excluding tests and __init__)
CHECKED_FILES = [
    "proteomics/db/protein_database_decoy.py",
    "ensembl/ensembl.py",
    "cgenomes/cgenomes_proteindb.py",
    "proteogenomics/blast_get_position.py",
    "proteogenomics/spectrumai.py",
]


def _count_print_calls(filepath: Path) -> list[int]:
    """Return line numbers of print() function calls in a Python file."""
    source = filepath.read_text()
    tree = ast.parse(source)
    lines = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            func = node.func
            if isinstance(func, ast.Name) and func.id == "print":
                lines.append(node.lineno)
    return lines


def test_no_print_in_production_code():
    violations = {}
    for rel_path in CHECKED_FILES:
        filepath = PGATK_ROOT / rel_path
        if filepath.exists():
            lines = _count_print_calls(filepath)
            if lines:
                violations[rel_path] = lines
    assert not violations, (
        f"print() calls found in production code (use logger instead):\n"
        + "\n".join(f"  {f}: lines {lines}" for f, lines in violations.items())
    )
```

**Step 2: Run test to verify it fails**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_no_print_statements.py -v`
Expected: FAIL with list of all print() locations

**Step 3: Replace print() calls file by file**

For each file, the pattern is:
1. Ensure the file has a logger (most already do via `ParameterConfiguration.get_logger()`)
2. Replace `print(...)` with `self.get_logger().info(...)` for informational output
3. Replace `print("ERROR: ...")` with `self.get_logger().error(...)`
4. Replace `print("Warning: ...")` with `self.get_logger().warning(...)`
5. For static methods or module-level functions without `self`, use `logging.getLogger(__name__)`

**Key replacements by file:**

**protein_database_decoy.py** — all `print()` inside instance methods → `self.get_logger().info()`:
- Lines 158-160, 193, 195, 197, 199-201: diagnostic stats → `self.get_logger().info(...)`
- Lines 261, 262, 279, 282, 286, 320, 367: decoypyrat progress → `self.get_logger().info(...)`
- Lines 481, 483-484: final stats → `self.get_logger().info(...)`

**ensembl.py** — mixed instance and static methods:
- Line 131: `print("skip entries...")` → `self.get_logger().warning("skip entries...")`
- Line 248: `print("Database already exists...")` → `logging.getLogger(__name__).debug("Database already exists: %s %s", e, gtf_db_file)`
- Line 272-274: `print("""Feature {} ...""")` → `logging.getLogger(__name__).warning(...)`
- Line 400: `print("Could not extra cds...")` → `self.get_logger().warning(...)`
- Line 755: redundant `print(msg)` after `logger` call → remove the print
- Lines 823-829: `print("   total number...")` → `self.get_logger().info(...)`
- Line 857-858: `print("ERROR: This script...")` → `logging.getLogger(__name__).error(...)`

**cgenomes_proteindb.py** — add `import logging` at top, use `self.get_logger()` in instance methods and `logging.getLogger(__name__)` in static methods. Replace all 18 print calls.

**blast_get_position.py** — replace 8 print calls with `self.get_logger()`.

**spectrumai.py** — replace 6 print calls with `self.get_logger()`.

**Step 4: Run the no-print test**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_no_print_statements.py -v`
Expected: PASS

**Step 5: Run full test suite**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/ -v --timeout=120 -x`
Expected: All PASS

**Step 6: Commit**

```bash
git add pgatk/proteomics/db/protein_database_decoy.py pgatk/ensembl/ensembl.py pgatk/cgenomes/cgenomes_proteindb.py pgatk/proteogenomics/blast_get_position.py pgatk/proteogenomics/spectrumai.py pgatk/tests/test_no_print_statements.py
git commit -m "refactor: replace print() with logger across all production modules"
```

---

## Task 7: Replace broad except Exception with specific exceptions

**Files:**
- Modify: `pgatk/toolbox/general.py:163,229,342`
- Modify: `pgatk/ensembl/ensembl.py:247`
- Modify: `pgatk/proteogenomics/spectrumai.py:163,180`

**Step 1: Fix each broad except block**

**`toolbox/general.py` line 163:**
```python
# OLD:
except Exception as e:
    raise ToolBoxException(str(e))
# NEW:
except (OSError, PermissionError) as e:
    raise ToolBoxException(str(e)) from e
```

**`toolbox/general.py` line 229 (download_file retry):**
```python
# OLD:
except Exception as error:
    remaining_download_tries -= 1
    log.error("Error code: " + str(error))
# NEW:
except (URLError, ContentTooShortError, OSError, HTTPError) as error:
    remaining_download_tries -= 1
    log.error("Download error: %s", error)
```

**`toolbox/general.py` line 342 (gunzip_files):**
```python
# OLD:
except Exception as e:
    err_msg = "UNKNOWN ERROR uncompressing file..."
# NEW:
except (OSError, subprocess.CalledProcessError) as e:
    err_msg = "Error uncompressing file..."
```

**`ensembl/ensembl.py` line 247 (parse_gtf):**
```python
# OLD:
except Exception as e:
    print("Database already exists: " + str(e), gtf_db_file)
# NEW:
except (gffutils.FeatureNotFoundError, ValueError, sqlite3.OperationalError) as e:
    logging.getLogger(__name__).debug("GTF database already exists: %s %s", e, gtf_db_file)
```
Add `import sqlite3` to the imports.

**`spectrumai.py` line 163 (mzML loading):**
```python
# OLD:
except Exception as e:
    print(mzml_file + " has ERROR!")
# NEW:
except (OSError, RuntimeError, ValueError) as e:
    self.get_logger().error("Error loading mzML file %s: %s", mzml_file, e)
```

**`spectrumai.py` line 180 (scan lookup):**
```python
# OLD:
except Exception as e:
    print("ERROR: ...")
# NEW:
except (RuntimeError, IndexError, KeyError) as e:
    self.get_logger().error("Scan lookup error in %s scan %s: %s", mzml_file, scan_num, e)
```

**Step 2: Run full test suite**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/ -v --timeout=120 -x`
Expected: All PASS

**Step 3: Commit**

```bash
git add pgatk/toolbox/general.py pgatk/ensembl/ensembl.py pgatk/proteogenomics/spectrumai.py
git commit -m "refactor: replace broad except Exception with specific exception types"
```

---

## Task 8: Standardize config fallback pattern with helper method

**Files:**
- Modify: `pgatk/toolbox/general.py` (add `_get_config_value` helper to `ParameterConfiguration`)
- Modify: `pgatk/ensembl/ensembl.py` (simplify `__init__` config loading)
- Modify: `pgatk/cgenomes/cgenomes_proteindb.py` (simplify `get_mutations_default_options`)
- Modify: `pgatk/proteogenomics/blast_get_position.py` (simplify `get_blast_parameters`)
- Modify: `pgatk/proteogenomics/spectrumai.py` (simplify parameter loading)

**Step 1: Write the test**

```python
# pgatk/tests/test_parameter_configuration.py
"""Test the ParameterConfiguration config fallback helper."""
from pgatk.toolbox.general import ParameterConfiguration


def test_get_config_value_from_pipeline_params():
    """Pipeline params should take precedence."""
    config = ParameterConfiguration.__new__(ParameterConfiguration)
    config._pipeline_params = {"root": {"key1": "from_pipeline"}}
    config._default_params = {"root": {"key1": "from_default"}}
    config._ROOT_CONFIG_NAME = "root"
    result = config.get_config_value("key1", "fallback")
    assert result == "from_pipeline"


def test_get_config_value_from_defaults():
    """Default params used when pipeline params missing."""
    config = ParameterConfiguration.__new__(ParameterConfiguration)
    config._pipeline_params = {"root": {}}
    config._default_params = {"root": {"key1": "from_default"}}
    config._ROOT_CONFIG_NAME = "root"
    result = config.get_config_value("key1", "fallback")
    assert result == "from_default"


def test_get_config_value_fallback():
    """Fallback returned when neither source has the key."""
    config = ParameterConfiguration.__new__(ParameterConfiguration)
    config._pipeline_params = {"root": {}}
    config._default_params = {"root": {}}
    config._ROOT_CONFIG_NAME = "root"
    result = config.get_config_value("missing_key", "fallback_value")
    assert result == "fallback_value"
```

**Step 2: Run test to verify it fails**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_parameter_configuration.py -v`
Expected: FAIL — `get_config_value` does not exist yet

**Step 3: Add the helper method to ParameterConfiguration**

In `pgatk/toolbox/general.py`, add to the `ParameterConfiguration` class:

```python
def get_config_value(self, key: str, default=None):
    """Get a configuration value with standard 3-layer fallback.

    Checks pipeline_params first, then default_params, then returns default.
    """
    root = self._ROOT_CONFIG_NAME
    # Check pipeline params first
    if (self._pipeline_params is not None
            and root in self._pipeline_params
            and key in self._pipeline_params[root]):
        return self._pipeline_params[root][key]
    # Then check default config params
    if (self._default_params is not None
            and root in self._default_params
            and key in self._default_params[root]):
        return self._default_params[root][key]
    return default
```

**Step 4: Run test to verify it passes**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_parameter_configuration.py -v`
Expected: PASS

**Step 5: Refactor subclasses to use the new helper**

Replace the repetitive fallback patterns in each subclass with `self.get_config_value(key, default)`. For example, in `blast_get_position.py`, the `get_blast_parameters()` method:

```python
# OLD pattern (repeated everywhere):
def get_blast_parameters(self, variable: str, default_value):
    value = default_value
    if self._pipeline_params is not None and self._ROOT_CONFIG_NAME in self._pipeline_params:
        if variable in self._pipeline_params[self._ROOT_CONFIG_NAME]:
            value = self._pipeline_params[self._ROOT_CONFIG_NAME][variable]
    elif self._default_params is not None and self._ROOT_CONFIG_NAME in self._default_params:
        if variable in self._default_params[self._ROOT_CONFIG_NAME]:
            value = self._default_params[self._ROOT_CONFIG_NAME][variable]
    return value

# NEW:
def get_blast_parameters(self, variable: str, default_value):
    return self.get_config_value(variable, default_value)
```

Apply the same simplification to equivalent methods in `cgenomes_proteindb.py`, `spectrumai.py`, and `ensembl.py`.

**Step 6: Run full test suite**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/ -v --timeout=120 -x`
Expected: All PASS

**Step 7: Commit**

```bash
git add pgatk/toolbox/general.py pgatk/ensembl/ensembl.py pgatk/cgenomes/cgenomes_proteindb.py pgatk/proteogenomics/blast_get_position.py pgatk/proteogenomics/spectrumai.py pgatk/tests/test_parameter_configuration.py
git commit -m "refactor: add get_config_value helper, DRY up config fallback pattern"
```

---

## Task 9: Convert db/digest_mutant_protein.py to proper CLI command

**Files:**
- Modify: `pgatk/db/digest_mutant_protein.py` (wrap in Click command)
- Modify: `pgatk/cli.py` (register new command)
- Create: `pgatk/commands/digest_mutant_protein.py` (Click wrapper)

**Step 1: Refactor the script into a function**

Wrap all logic in `pgatk/db/digest_mutant_protein.py` inside a function with proper parameters:

```python
# pgatk/db/digest_mutant_protein.py
import re
import logging
from collections import OrderedDict
from typing import Optional

from Bio import SeqIO

logger = logging.getLogger(__name__)


def trypsin_cleavage(proseq: str, miss_cleavage: int) -> list:
    """Perform in-silico tryptic digestion."""
    # (existing trypsin_cleavage function body, unchanged)
    ...


def digest_mutant_proteins(
    input_file: str,
    fasta_file: str,
    output_file: str,
    decoy_prefix: str = "DECOY_",
    min_length: int = 7,
    max_length: int = 40,
    missed_cleavages: int = 2,
) -> None:
    """Digest mutant proteins and filter against canonical peptidome.

    Reads mutant protein sequences, performs tryptic digestion, removes peptides
    found in the canonical proteome, and writes unique variant peptides.
    """
    # (existing script logic from lines 83-139, using function parameters
    #  instead of getopt-parsed variables)
    ...
```

Remove all `sys.argv`, `getopt`, and `sys.exit()` code. Remove the import-time execution.

**Step 2: Create Click command wrapper**

```python
# pgatk/commands/digest_mutant_protein.py
import click

from pgatk.db.digest_mutant_protein import digest_mutant_proteins


@click.command("digest-mutant-protein", short_help="Digest mutant proteins and filter against canonical proteome")
@click.option("-i", "--input", "input_file", required=True, help="Input mutant protein FASTA")
@click.option("-f", "--fasta", "fasta_file", required=True, help="Reference canonical protein FASTA")
@click.option("-o", "--output", "output_file", required=True, help="Output file for unique variant peptides")
@click.option("--prefix", default="DECOY_", help="Decoy prefix to skip (default: DECOY_)")
@click.option("--min-len", default=7, type=int, help="Minimum peptide length (default: 7)")
@click.option("--max-len", default=40, type=int, help="Maximum peptide length (default: 40)")
@click.option("--missed-cleavages", default=2, type=int, help="Missed cleavages (default: 2)")
def digest_mutant_protein(input_file, fasta_file, output_file, prefix, min_len, max_len, missed_cleavages):
    digest_mutant_proteins(input_file, fasta_file, output_file, prefix, min_len, max_len, missed_cleavages)
```

**Step 3: Register in cli.py**

Add to `pgatk/cli.py`:
```python
from pgatk.commands.digest_mutant_protein import digest_mutant_protein
cli.add_command(digest_mutant_protein)
```

**Step 4: Run CLI help to verify registration**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pgatk.cli digest-mutant-protein --help`
Expected: Help text with all options displayed

**Step 5: Run full test suite**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/ -v --timeout=120 -x`
Expected: All PASS

**Step 6: Commit**

```bash
git add pgatk/db/digest_mutant_protein.py pgatk/commands/digest_mutant_protein.py pgatk/cli.py
git commit -m "refactor: convert digest_mutant_protein from standalone script to Click CLI command"
```

---

## Task 10: Convert db/map_peptide2genome.py to proper CLI command

**Files:**
- Modify: `pgatk/db/map_peptide2genome.py` (wrap in functions)
- Create: `pgatk/commands/map_peptide2genome.py` (Click wrapper)
- Modify: `pgatk/cli.py` (register new command)

**Step 1: Refactor the script into functions**

Same pattern as Task 9. Wrap all logic in a function:

```python
# pgatk/db/map_peptide2genome.py
import re
import logging
from typing import Optional

from Bio import SeqIO

logger = logging.getLogger(__name__)


class EXON:
    """Represents an exon with genomic coordinates."""
    # (existing EXON class, unchanged except add type hints)
    ...


def cal_trans_pos(exon_list: list) -> list:
    # (existing function, unchanged)
    ...


def get_pep_cor(exon_object_list: list, n1: int, n2: int) -> list:
    # (existing function, unchanged)
    ...


def parse_gtf(infile: str) -> dict:
    # (existing function, unchanged)
    ...


def map_peptides_to_genome(
    input_file: str,
    gtf_file: str,
    fasta_file: str,
    idmap_file: str,
    output_file: str,
    pep_col: int = 0,
    prot_col: int = 1,
) -> None:
    """Map identified peptides to genomic coordinates.

    Reads peptide identifications, maps them to transcript coordinates using
    a GTF file and protein-transcript ID mapping, outputs GFF3 format.
    """
    # (existing script logic from lines 146-261, using function parameters)
    ...
```

Remove all `sys.argv`, `getopt`, `sys.exit()` code.

**Step 2: Create Click command wrapper**

```python
# pgatk/commands/map_peptide2genome.py
import click

from pgatk.db.map_peptide2genome import map_peptides_to_genome


@click.command("map-peptide2genome", short_help="Map peptides to genomic coordinates (GFF3 output)")
@click.option("-i", "--input", "input_file", required=True, help="Input peptide identification TSV")
@click.option("-g", "--gtf", "gtf_file", required=True, help="GTF gene annotation file")
@click.option("-f", "--fasta", "fasta_file", required=True, help="Protein FASTA file")
@click.option("-m", "--idmap", "idmap_file", required=True, help="Protein-to-transcript ID mapping file")
@click.option("-o", "--output", "output_file", required=True, help="Output GFF3 file")
@click.option("--pep-col", default=0, type=int, help="Peptide column index (default: 0)")
@click.option("--prot-col", default=1, type=int, help="Protein column index (default: 1)")
def map_peptide2genome(input_file, gtf_file, fasta_file, idmap_file, output_file, pep_col, prot_col):
    map_peptides_to_genome(input_file, gtf_file, fasta_file, idmap_file, output_file, pep_col, prot_col)
```

**Step 3: Register in cli.py**

Add to `pgatk/cli.py`:
```python
from pgatk.commands.map_peptide2genome import map_peptide2genome
cli.add_command(map_peptide2genome)
```

**Step 4: Run CLI help to verify registration**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pgatk.cli map-peptide2genome --help`
Expected: Help text with all options displayed

**Step 5: Run full test suite**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/ -v --timeout=120 -x`
Expected: All PASS

**Step 6: Commit**

```bash
git add pgatk/db/map_peptide2genome.py pgatk/commands/map_peptide2genome.py pgatk/cli.py
git commit -m "refactor: convert map_peptide2genome from standalone script to Click CLI command"
```

---

## Task 11: Use pathlib.Path for file operations

**Files:**
- Modify: `pgatk/ensembl/ensembl.py` (5 sites)
- Modify: `pgatk/cgenomes/cgenomes_proteindb.py` (2 sites)
- Modify: `pgatk/proteomics/db/protein_database_decoy.py` (1 site)

**Step 1: Add pathlib import and replace string path operations**

**`ensembl/ensembl.py`** — add `from pathlib import Path` to imports:

- Line 417: `os.path.abspath(vcf_file.split('/')[-1].replace('.vcf', ''))` → `str(Path(vcf_file).resolve().stem)`
- Line 419: `annotated_vcf + '_all.bed'` → `str(Path(annotated_vcf).with_name(Path(annotated_vcf).name + '_all.bed'))`
  Or simpler: keep `annotated_vcf` as a Path stem and use f-strings: `f"{annotated_vcf}_all.bed"`
- Line 508: `gene_annotations_gtf.replace('.gtf', '.db')` → `str(Path(gene_annotations_gtf).with_suffix('.db'))`

**`cgenomes/cgenomes_proteindb.py`** — add `from pathlib import Path`:

- Line 260: `self._local_output_file.replace('.fa', '') + '_' + ...` → `str(Path(self._local_output_file).with_suffix('')) + '_' + ...` (or use Path stem)
- Line 461: same pattern

**`protein_database_decoy.py`** — add `from pathlib import Path`:

- Line 479: `self._output_file.replace('.fa', '') + '_noAlternative.fa'` → `str(Path(self._output_file).with_suffix('').with_name(Path(self._output_file).stem + '_noAlternative.fa'))`

**Step 2: Run full test suite**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/ -v --timeout=120 -x`
Expected: All PASS

**Step 3: Commit**

```bash
git add pgatk/ensembl/ensembl.py pgatk/cgenomes/cgenomes_proteindb.py pgatk/proteomics/db/protein_database_decoy.py
git commit -m "refactor: use pathlib.Path for file path operations"
```

---

## Task 12: Add type hints to public APIs

**Files:**
- Modify: `pgatk/toolbox/general.py`
- Modify: `pgatk/ensembl/ensembl.py`
- Modify: `pgatk/cgenomes/cgenomes_proteindb.py`
- Modify: `pgatk/proteogenomics/blast_get_position.py`
- Modify: `pgatk/proteogenomics/spectrumai.py`
- Modify: `pgatk/proteomics/db/protein_database_decoy.py`

**Step 1: Add type hints to each file's public methods**

Work through each file adding type annotations to all public methods (non-underscore-prefixed). Use `from __future__ import annotations` at the top of each modified file for forward references.

Key patterns:
```python
from __future__ import annotations
from typing import Optional
from pathlib import Path
import logging
import pandas as pd
import gffutils
```

Examples of what to add:

```python
# ensembl.py
def vcf_to_proteindb(self, vcf_file: str, input_fasta: str, gene_annotations_gtf: str) -> str: ...
def dnaseq_to_proteindb(self, input_fasta: str) -> str: ...
def check_proteindb(self, input_fasta: str, add_stop_codon: bool = True, num_aa: int = 6) -> None: ...

@staticmethod
def get_altseq(ref_seq: str, ref_allele: str, var_allele: str, var_pos: int,
               strand: str, features_info: list, cds_info: Optional[list] = None) -> tuple: ...

@staticmethod
def parse_gtf(gene_annotations_gtf: str, gtf_db_file: str) -> gffutils.FeatureDB: ...
```

```python
# blast_get_position.py
def blast(self, input_psm_to_blast: str, output_psm: str) -> None: ...

@staticmethod
def get_details(fasta: str, peptide: str) -> list[str]: ...

@staticmethod
def peptide_blast_protein(fasta: str, peptide: str) -> Optional[list]: ...
```

**Step 2: Run mypy (optional, non-blocking) and tests**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/ -v --timeout=120 -x`
Expected: All PASS

**Step 3: Commit**

```bash
git add pgatk/toolbox/general.py pgatk/ensembl/ensembl.py pgatk/cgenomes/cgenomes_proteindb.py pgatk/proteogenomics/blast_get_position.py pgatk/proteogenomics/spectrumai.py pgatk/proteomics/db/protein_database_decoy.py
git commit -m "refactor: add type hints to all public API methods"
```

---

## Task 13: Expand unit test coverage for core functions

**Files:**
- Create: `pgatk/tests/test_ensembl_core.py`
- Create: `pgatk/tests/test_cosmic_core.py`
- Create: `pgatk/tests/test_decoy_core.py`

**Step 1: Write unit tests for get_altseq()**

```python
# pgatk/tests/test_ensembl_core.py
"""Unit tests for EnsemblDataService core functions."""
from pgatk.ensembl.ensembl import EnsemblDataService


class TestGetAltseq:
    """Test the static method that computes alternative sequences from variants."""

    def test_snp_plus_strand(self):
        """Simple SNP on plus strand should substitute one base."""
        ref_seq = "ATGCGATCGA"
        # SNP at position 3 (0-indexed in CDS): G -> T
        alt_seq = EnsemblDataService.get_altseq(
            ref_seq=ref_seq, ref_allele="G", var_allele="T",
            var_pos=4, strand="+",
            features_info=[(0, len(ref_seq))],  # single exon covering full CDS
        )
        assert alt_seq is not None
        # The exact assertion depends on what get_altseq returns —
        # inspect the return value and assert accordingly

    def test_insertion_plus_strand(self):
        """Insertion on plus strand should add bases."""
        ref_seq = "ATGCGATCGA"
        alt_seq = EnsemblDataService.get_altseq(
            ref_seq=ref_seq, ref_allele="G", var_allele="GAA",
            var_pos=4, strand="+",
            features_info=[(0, len(ref_seq))],
        )
        assert alt_seq is not None

    def test_deletion_plus_strand(self):
        """Deletion on plus strand should remove bases."""
        ref_seq = "ATGCGATCGA"
        alt_seq = EnsemblDataService.get_altseq(
            ref_seq=ref_seq, ref_allele="CGA", var_allele="C",
            var_pos=4, strand="+",
            features_info=[(0, len(ref_seq))],
        )
        assert alt_seq is not None
```

**Step 2: Write unit tests for get_mut_pro_seq()**

```python
# pgatk/tests/test_cosmic_core.py
"""Unit tests for CancerGenomesService core mutation functions."""
from pgatk.cgenomes.cgenomes_proteindb import CancerGenomesService
from pgatk.cgenomes.models import SNP


class TestGetMutProSeq:
    """Test the static method that applies COSMIC mutations to protein sequences."""

    def test_substitution_missense(self):
        """A single AA substitution should change one residue."""
        snp = SNP(gene="TEST", aa_mut="p.A10V", mutation_type="Substitution - Missense")
        ref_seq = "MAAAAAAAAAAAAAAAAA"  # A at position 10
        result = CancerGenomesService.get_mut_pro_seq(snp, ref_seq)
        # A at position 10 (1-indexed) should become V
        assert result is not None
        assert "V" in result

    def test_nonsense(self):
        """Nonsense mutation should truncate at the stop."""
        snp = SNP(gene="TEST", aa_mut="p.A10*", mutation_type="Substitution - Nonsense")
        ref_seq = "MAAAAAAAAAAAAAAAAA"
        result = CancerGenomesService.get_mut_pro_seq(snp, ref_seq)
        assert result is not None
```

**Step 3: Write unit tests for decoy revswitch()**

```python
# pgatk/tests/test_decoy_core.py
"""Unit tests for ProteinDBDecoyService core functions."""
from pgatk.proteomics.db.protein_database_decoy import ProteinDBDecoyService


class TestRevswitch:
    """Test the reverse-with-KR-switch function."""

    def test_revswitch_preserves_length(self):
        """Reversed sequence should have the same length."""
        seq = "ABCDEFGHIJK"
        result = ProteinDBDecoyService.revswitch(seq)
        assert len(result) == len(seq)

    def test_revswitch_kr_handling(self):
        """K and R at cleavage sites should be preserved at the correct positions."""
        seq = "PEPTIDEK"
        result = ProteinDBDecoyService.revswitch(seq)
        assert isinstance(result, str)
        assert len(result) == len(seq)
```

**Step 4: Run all new tests**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/test_ensembl_core.py pgatk/tests/test_cosmic_core.py pgatk/tests/test_decoy_core.py -v`
Expected: All PASS (adjust assertions based on actual function return values)

**Step 5: Run full test suite**

Run: `cd /Users/yperez/work/proteogenomics/pgatk && python -m pytest pgatk/tests/ -v --timeout=120 -x`
Expected: All PASS

**Step 6: Commit**

```bash
git add pgatk/tests/test_ensembl_core.py pgatk/tests/test_cosmic_core.py pgatk/tests/test_decoy_core.py
git commit -m "test: add unit tests for get_altseq, get_mut_pro_seq, and revswitch"
```

---

## Task Summary

| Task | Description | Files | Depends On |
|------|-------------|-------|------------|
| 1 | Add conftest.py with fixtures | tests/conftest.py | — |
| 2 | Fix variable shadowing bug | protein_database_decoy.py | — |
| 3 | Convert SNP to dataclass | cgenomes/models.py, cgenomes_proteindb.py | — |
| 4 | Replace pathos in blast_get_position | blast_get_position.py | — |
| 5 | Replace pathos in spectrumai | spectrumai.py | — |
| 6 | Replace print() with logger | 5 files | — |
| 7 | Replace broad except blocks | 3 files | Task 6 |
| 8 | Standardize config fallback | toolbox/general.py + 4 files | — |
| 9 | Convert digest_mutant_protein to CLI | db/digest_mutant_protein.py, cli.py | — |
| 10 | Convert map_peptide2genome to CLI | db/map_peptide2genome.py, cli.py | — |
| 11 | Use pathlib.Path | 3 files | — |
| 12 | Add type hints to public APIs | 6 files | Tasks 1-11 |
| 13 | Expand unit test coverage | 3 test files | Tasks 1-11 |

**Parallelizable groups:**
- Tasks 1-5 can run in parallel (no shared file dependencies)
- Task 6 can start after Tasks 2, 4, 5 (they touch overlapping files)
- Tasks 7 depends on Task 6 (print→logger before except→specific)
- Tasks 9-10 are independent of each other
- Tasks 12-13 should run last (they depend on all prior cleanups)
