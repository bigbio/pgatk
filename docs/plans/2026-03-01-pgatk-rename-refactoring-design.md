# Design: py-pgatk → pgatk Rename & Refactoring

**Date:** 2026-03-01
**Status:** Approved

## Context

The project is currently named `py-pgatk` with a Python package `pypgatk`. We are doing a full cleanup: rename, config system overhaul, bug fixes, and dead code removal.

## Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Scope | Full cleanup | Rename + config + bugs + dead code in one effort |
| Delivery | 4 separate PRs | Reviewable, independent after PR 1 |
| Backward compat | Clean break | No shim package, rename everything |
| Config approach | Config Registry | Central mapping, auto-load bundled defaults |
| CLI entry point | `pgatk` | Matches project name |

## PR 1: Rename (branch: `refactor/rename`)

**Rename `pypgatk/` → `pgatk/`** across the entire codebase.

### Changes

- Rename directory: `pypgatk/` → `pgatk/`
- Rename CLI module: `pypgatkc.py` → `cli.py`
- Update `pyproject.toml`:
  - `name = "pgatk"`
  - Entry point: `pgatk = "pgatk.cli:main"`
  - Package find: `include = ["pgatk*"]`
  - Package data: `pgatk = ["config/*.yaml", "config/*.json"]`
  - URLs: `bigbio/pgatk`
- Update all internal imports: `from pypgatk.xxx` → `from pgatk.xxx`
- Update docs references (mkdocs.yml, README.md, doc pages)
- Update GitHub Actions workflows
- Update test imports

### Files affected

- Every `.py` file with `pypgatk` imports (~44 files)
- `pyproject.toml`
- `mkdocs.yml`
- `README.md`
- `.github/workflows/*.yml`
- `docs/*.md`

## PR 2: Config Registry (branch: `refactor/config-registry`)

**Auto-load bundled config defaults so `-c` flag becomes optional.**

### New module: `pgatk/config/registry.py`

```python
import os
from pgatk.toolbox.general import read_yaml_from_file

CONFIG_DIR = os.path.dirname(__file__)

COMMAND_CONFIGS = {
    "ensembl_downloader":    "ensembl_downloader_config.yaml",
    "ensembl_config":        "ensembl_config.yaml",
    "cosmic":                "cosmic_config.yaml",
    "cbioportal":            "cbioportal_config.yaml",
    "protein_decoy":         "protein_decoy.yaml",
    "openms_analysis":       "openms_analysis.yaml",
}

def load_config(config_key, user_config_file=None):
    """Load config: user file if provided, otherwise bundled default."""
    if user_config_file is not None:
        return read_yaml_from_file(user_config_file)
    default_file = COMMAND_CONFIGS.get(config_key)
    if default_file is None:
        return {}
    return read_yaml_from_file(os.path.join(CONFIG_DIR, default_file))
```

### Command changes

Replace the 3-line boilerplate in every command:

```python
# Before
config_data = None
if config_file is not None:
    config_data = read_yaml_from_file(config_file)

# After
from pgatk.config.registry import load_config
config_data = load_config("ensembl_downloader", config_file)
```

### Files affected

- New: `pgatk/config/registry.py`
- All 14 command files in `pgatk/commands/`

## PR 3: Bug Fixes (branch: `fix/bugs`)

| Bug | File | Fix |
|---|---|---|
| `_filter_write_idxml_with_df` writes unfiltered data | `proteomics/openms.py` | Write only filtered entries |
| ISO-8859 encoding typo | `cgenomes/cosmic_downloader.py` | `'ISO-8859'` → `'ISO-8859-1'` |
| Deprecated `Bio.pairwise2` | `proteomics/openms.py` | Replace with `Bio.Align.PairwiseAligner` |
| Enzyme name typos | `proteogenomics/models.py` | `'ArgC/P'` → `'Arg-C/P'` etc. |
| `print_help()` for required args | `commands/*.py` | Use Click's `required=True` |
| `download_file` log overwrite | `toolbox/general.py:193` | `if log is not None:` → `if log is None:` |

## PR 4: Cleanup (branch: `refactor/cleanup`)

| Item | File | Action |
|---|---|---|
| Dead `Species`/`SpeciesService` classes | `ensembl/ensembl.py` | Remove if unused |
| Commented-out code after return | `toolbox/general.py:241-253` | Delete |
| Hardcoded version `'0.0.26'` | `cli.py` | Use `importlib.metadata.version("pgatk")` |
| Log files created in CWD | `toolbox/general.py` | Log to stderr only by default |

## Dependency Order

```
PR 1 (Rename) ──► PR 2 (Config Registry)
                ──► PR 3 (Bug Fixes)
                ──► PR 4 (Cleanup)
```

PR 1 must merge first. PRs 2, 3, 4 are independent of each other.
