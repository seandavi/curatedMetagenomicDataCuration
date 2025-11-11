# UV Quick Reference Guide

This guide shows how to use [uv](https://github.com/astral-sh/uv) with the Metadata Validator Framework.

## What is uv?

uv is an extremely fast Python package installer and resolver, written in Rust. It's 10-100x faster than pip and has better dependency resolution.

## Installation

### Install uv

```bash
# On macOS/Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# On Windows
powershell -c "irm https://astral.sh/uv/install.ps1 | iex"
```

## Common Workflows

### 1. Initial Setup

```bash
# Navigate to the project
cd metadata_validator

# Create a virtual environment
uv venv

# Activate the virtual environment
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install the package in editable mode
uv pip install -e .
```

### 2. Install with Optional Dependencies

```bash
# Install with ontology support (pronto, requests)
uv pip install -e ".[ontology]"

# Install with fuzzy matching support
uv pip install -e ".[fuzzy]"

# Install with development tools
uv pip install -e ".[dev]"

# Install everything
uv pip install -e ".[all]"
```

### 3. Running Scripts

```bash
# After activating venv
python test_validator.py
python examples/simple_example.py

# Or use uv run (no need to activate venv)
uv run python test_validator.py
uv run python examples/simple_example.py
```

### 4. Adding New Dependencies

If you need to add new packages:

```bash
# Add to pyproject.toml dependencies, then:
uv pip install -e .

# Or install directly
uv pip install package-name
```

### 5. Updating Dependencies

```bash
# Update all packages
uv pip install --upgrade -e .

# Update a specific package
uv pip install --upgrade package-name
```

### 6. List Installed Packages

```bash
uv pip list
```

### 7. Remove Virtual Environment

```bash
# Deactivate if active
deactivate

# Remove the directory
rm -rf .venv
```

## Key Advantages of uv

1. **Speed**: 10-100x faster than pip
2. **Better Dependency Resolution**: More reliable conflict resolution
3. **Reproducible**: Consistent installs across machines
4. **Modern**: Built with Rust, actively developed
5. **Compatible**: Drop-in replacement for pip

## Common Commands Comparison

| Task | pip | uv |
|------|-----|-----|
| Create venv | `python -m venv .venv` | `uv venv` |
| Install package | `pip install package` | `uv pip install package` |
| Install from pyproject.toml | `pip install -e .` | `uv pip install -e .` |
| Install extras | `pip install -e ".[dev]"` | `uv pip install -e ".[dev]"` |
| List packages | `pip list` | `uv pip list` |
| Freeze deps | `pip freeze` | `uv pip freeze` |

## Troubleshooting

### uv command not found

After installation, restart your terminal or source your shell config:

```bash
source ~/.bashrc  # or ~/.zshrc
```

### Python version issues

uv uses the `.python-version` file in this directory (currently set to 3.11). To use a different version:

```bash
echo "3.10" > .python-version
uv venv
```

### Import errors

Make sure you've activated the virtual environment or use `uv run`:

```bash
source .venv/bin/activate
python test_validator.py
```

## CI/CD Integration

Example GitHub Actions workflow:

```yaml
name: Test
on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3

      - name: Install uv
        run: curl -LsSf https://astral.sh/uv/install.sh | sh

      - name: Install dependencies
        working-directory: metadata_validator
        run: |
          uv venv
          source .venv/bin/activate
          uv pip install -e ".[dev]"

      - name: Run tests
        working-directory: metadata_validator
        run: |
          source .venv/bin/activate
          python test_validator.py
```

## Learn More

- [uv Documentation](https://github.com/astral-sh/uv)
- [uv Installation Guide](https://github.com/astral-sh/uv#installation)
- [Astral (uv creators)](https://astral.sh/)

## Back to pip?

If you prefer to use pip, no problem! The `pyproject.toml` works with both:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```
