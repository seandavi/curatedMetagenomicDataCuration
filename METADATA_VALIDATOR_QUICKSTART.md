# Metadata Validator Framework - Quick Start

## Overview

The Metadata Validator Framework provides comprehensive validation and harmonization capabilities for sample metadata using an agentic architecture. The framework is located in the `metadata_validator/` directory.

## Installation

### Using uv (Recommended - Fast!)

```bash
# Install uv if you haven't already
curl -LsSf https://astral.sh/uv/install.sh | sh

# Navigate to the validator directory
cd metadata_validator

# Create a virtual environment and install
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
uv pip install -e .

# Or install directly without creating a venv
uv pip install -e .
```

### Using pip (Alternative)

```bash
cd metadata_validator
pip install -e .
```

## Quick Start Example

```python
from metadata_validator import MetadataOrchestrator

# Create orchestrator
orchestrator = MetadataOrchestrator(
    data_dictionary_path="inst/extdata/cMD_data_dictionary.csv",
    ontology_backend='custom',
    human_review=True,
    auto_approve=False
)

# Process a metadata file
summary = orchestrator.process_file(
    input_csv_path="your_metadata.csv",
    output_csv_path="your_metadata_harmonized.csv",
    report_path="validation_report.txt",
    interactive_review=True
)

# View summary
print(orchestrator.generate_summary_report(summary))
```

## Testing

Run the test suite to verify installation:

```bash
cd metadata_validator

# With uv (if using venv)
source .venv/bin/activate
python test_validator.py

# Or run directly with uv
uv run python test_validator.py
```

## Key Features

### 1. Validation
- Validates CSV metadata against the cMD data dictionary
- Supports regex patterns and enumerated values
- Handles multi-valued columns
- Generates detailed validation reports

### 2. Agentic Harmonization
Three specialized agents work together:

- **Column Renamer Agent**: Suggests column name corrections using fuzzy matching
- **Value Harmonizer Agent**: Fixes incorrect values to match allowed terms
- **Ontology Matcher Agent**: Maps terms to ontology concepts (NCIT, DOID, etc.)

### 3. Human-in-the-Loop
- Interactive review of proposed changes
- Batch approval strategies
- Complete provenance tracking

### 4. Workflow Orchestration
- End-to-end automated pipeline
- Keeps original values for audit trail
- Re-validates after harmonization
- Generates comprehensive reports

## Examples

See the `metadata_validator/examples/` directory for complete examples:

- `simple_example.py` - Basic workflow
- `advanced_example.py` - Custom configurations and advanced features

## Documentation

Full documentation is available in `metadata_validator/README.md`

## Architecture

```
metadata_validator/
├── validator.py              # Core validation engine
├── orchestrator.py           # Workflow coordinator
├── agents/                   # Specialized agents
│   ├── column_renamer_agent.py
│   ├── value_harmonizer_agent.py
│   ├── ontology_matcher_agent.py
│   └── human_review_agent.py
├── examples/                 # Example scripts
├── test_validator.py         # Test suite
├── pyproject.toml            # Project config & dependencies (uv/pip)
├── .python-version           # Python version for uv
└── README.md                 # Full documentation
```

## Output Files

The framework generates:

1. **Harmonized CSV**: Corrected metadata with original columns preserved
2. **Validation Report**: Detailed validation results
3. **Workflow History**: Complete record of all actions (JSON)
4. **Ontology Columns**: New columns with ontology term mappings

## Configuration

### Auto-approval Mode (for pipelines)
```python
orchestrator = MetadataOrchestrator(
    data_dictionary_path="inst/extdata/cMD_data_dictionary.csv",
    human_review=False,
    auto_approve=True
)
```

### Custom Ontology Backend
```python
# Option 1: Custom mappings (default)
orchestrator = MetadataOrchestrator(
    data_dictionary_path="...",
    ontology_backend='custom'
)

# Option 2: Local OBO file
orchestrator = MetadataOrchestrator(
    data_dictionary_path="...",
    ontology_backend='pronto',
    ontology_config={'ontology_path': '/path/to/ontology.obo'}
)

# Option 3: BioPortal API
orchestrator = MetadataOrchestrator(
    data_dictionary_path="...",
    ontology_backend='bioportal',
    ontology_config={
        'api_key': 'your-api-key',
        'ontologies': ['NCIT', 'DOID', 'HP']
    }
)
```

## Support

For issues or questions, please refer to the full documentation in `metadata_validator/README.md` or examine the example scripts in `metadata_validator/examples/`.
