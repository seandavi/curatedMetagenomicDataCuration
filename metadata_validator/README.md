# Metadata Validator Framework

A comprehensive Python framework for validating, harmonizing, and enriching sample metadata against a data dictionary using an agentic architecture.

## Features

### Core Validation
- **Data Dictionary Validation**: Validates CSV metadata files against a structured data dictionary
- **Regex Pattern Matching**: Supports both enumerated values and regular expression patterns
- **Multiple Value Support**: Handles columns with multiple delimited values
- **Comprehensive Reporting**: Generates detailed validation reports with error locations

### Agentic Harmonization System

The framework includes specialized agents that work together to improve metadata quality:

1. **Column Renamer Agent**: Suggests column name corrections using fuzzy string matching
2. **Value Harmonizer Agent**: Fixes incorrect values by matching to allowed terms
3. **Ontology Matcher Agent**: Maps metadata terms to ontology concepts (NCIT, DOID, etc.)
4. **Human Review Agent**: Provides interactive approval workflow for proposed changes

### Workflow Orchestration
- **Automated Pipeline**: End-to-end validation and harmonization workflow
- **Provenance Tracking**: Keeps original columns/values for audit trail
- **Action History**: Records all proposed and applied changes
- **Re-validation**: Automatically re-validates after harmonization

## Installation

### Using uv (Recommended)

[uv](https://github.com/astral-sh/uv) is a fast Python package installer and resolver. It's 10-100x faster than pip!

**📖 See [UV_GUIDE.md](UV_GUIDE.md) for a complete guide to using uv with this project.**

```bash
# Install uv if you haven't already
curl -LsSf https://astral.sh/uv/install.sh | sh

# Install the package and dependencies
cd metadata_validator
uv pip install -e .

# Or create a virtual environment and install
uv venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
uv pip install -e .
```

### Optional Dependencies

Install with ontology support:
```bash
uv pip install -e ".[ontology]"
```

Install with fuzzy matching:
```bash
uv pip install -e ".[fuzzy]"
```

Install all optional dependencies:
```bash
uv pip install -e ".[all]"
```

Install development dependencies:
```bash
uv pip install -e ".[dev]"
```

### Using pip (Alternative)

```bash
cd metadata_validator
pip install -e .

# With optional dependencies
pip install -e ".[all]"
```

## Quick Start

### Simple Usage with Orchestrator

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
    input_csv_path="sample_metadata.csv",
    output_csv_path="sample_metadata_harmonized.csv",
    report_path="validation_report.txt",
    interactive_review=True
)

# View summary
print(orchestrator.generate_summary_report(summary))
```

### Validation Only

```python
from metadata_validator import MetadataValidator

# Initialize validator
validator = MetadataValidator("inst/extdata/cMD_data_dictionary.csv")

# Validate a file
results = validator.validate_file("sample_metadata.csv")

# Generate report
report = validator.generate_report(results)
print(report)
```

## Data Dictionary Format

The data dictionary CSV must contain the following columns:

- **ColName**: Column name in the metadata
- **ColClass**: Data type (character, integer, double, numeric)
- **Unique**: Whether values must be unique (unique, non-unique)
- **Required**: Whether column is required (required, optional)
- **MultipleValues**: Whether multiple values allowed (TRUE, FALSE)
- **Description**: Human-readable description
- **AllowedValues**: Pipe-separated list of values OR regex pattern
- **Delimiter**: Delimiter for multiple values (e.g., ";")
- **Separater**: Additional separator information
- **DynamicEnum**: Dynamic enumeration flag
- **DynamicEnumProperty**: Dynamic enumeration property

### Example Data Dictionary Entry

```csv
ColName,ColClass,Unique,Required,MultipleValues,Description,AllowedValues,Delimiter
sex,character,non-unique,optional,FALSE,Biological sex,Female|Male|NA,NA
age,integer,non-unique,optional,FALSE,Age in years,[0-9]+,NA
disease,character,non-unique,optional,TRUE,Disease conditions,Healthy|Diabetes|...,;
```

## Agent System

### Column Renamer Agent

Suggests column name corrections for unrecognized columns:

```python
from metadata_validator import ColumnRenamerAgent, MetadataValidator

validator = MetadataValidator("data_dictionary.csv")
renamer = ColumnRenamerAgent(
    validator.data_dictionary,
    similarity_threshold=0.6  # Adjust sensitivity
)

# Analyze and suggest renamings
actions = renamer.analyze(df, validation_results)

# Review suggestions
for action in actions:
    print(action.description)
    print(f"  Similarity: {action.details['similarity_score']:.2f}")

# Apply approved actions
df_renamed = renamer.apply_actions(df, actions)
```

### Value Harmonizer Agent

Fixes incorrect values to match allowed terms:

```python
from metadata_validator import ValueHarmonizerAgent

harmonizer = ValueHarmonizerAgent(
    validator.data_dictionary,
    similarity_threshold=0.8
)

# Teach custom mappings
harmonizer.learn_mapping('sex', 'M', 'Male')
harmonizer.learn_mapping('sex', 'F', 'Female')

# Analyze and suggest corrections
actions = harmonizer.analyze(df, validation_results)

# Apply approved actions
df_harmonized = harmonizer.apply_actions(df, actions)

# Save learned mappings for reuse
harmonizer.export_learned_mappings("learned_mappings.json")
```

### Ontology Matcher Agent

Maps terms to ontology concepts:

```python
from metadata_validator import OntologyMatcherAgent

# Option 1: Custom backend (built-in mappings)
matcher = OntologyMatcherAgent(
    validator.data_dictionary,
    backend='custom'
)

# Option 2: Local OBO file
matcher = OntologyMatcherAgent(
    validator.data_dictionary,
    backend='pronto',
    ontology_config={'ontology_path': '/path/to/ontology.obo'}
)

# Option 3: BioPortal API
matcher = OntologyMatcherAgent(
    validator.data_dictionary,
    backend='bioportal',
    ontology_config={
        'api_key': 'your-api-key',
        'ontologies': ['NCIT', 'DOID', 'HP']
    }
)

# Add custom mapping
matcher.add_custom_mapping(
    term='Custom Disease',
    ontology='NCIT',
    term_id='NCIT:C12345',
    term_label='Custom Disease',
    definition='A custom disease definition'
)

# Analyze and add ontology columns
actions = matcher.analyze(df)
df_with_ontology = matcher.apply_actions(df, actions)
```

### Human Review Agent

Interactive review of proposed changes:

```python
from metadata_validator import HumanReviewAgent

reviewer = HumanReviewAgent(validator.data_dictionary)

# Interactive review
approved_actions = reviewer.review_actions(
    all_actions,
    interactive=True
)

# Batch approval strategies
reviewer.batch_approve_by_type(all_actions, 'add_ontology_mappings')
reviewer.batch_approve_by_agent(all_actions, 'ColumnRenamerAgent')

# Generate approval report
report = reviewer.generate_approval_report(all_actions)
print(report)
```

## Workflow Orchestration

The orchestrator coordinates all agents in a complete workflow:

### Workflow Steps

1. **Load Data**: Read input CSV file
2. **Initial Validation**: Validate against data dictionary
3. **Agent Analysis**: Each agent proposes corrections
4. **Human Review**: Interactive approval of proposed actions (optional)
5. **Apply Actions**: Execute approved changes in order:
   - Column renaming
   - Value harmonization
   - Ontology mapping
6. **Re-validation**: Validate harmonized data
7. **Save Results**: Export harmonized CSV and reports

### Workflow Configuration

```python
orchestrator = MetadataOrchestrator(
    data_dictionary_path="data_dictionary.csv",
    ontology_backend='custom',  # or 'pronto', 'bioportal'
    ontology_config={},  # Backend-specific config
    human_review=True,  # Enable human review
    auto_approve=False  # Require manual approval
)
```

### Auto-approval Mode

For automated pipelines:

```python
orchestrator = MetadataOrchestrator(
    data_dictionary_path="data_dictionary.csv",
    human_review=False,
    auto_approve=True
)

# All actions will be automatically approved and applied
summary = orchestrator.process_file(
    input_csv_path="input.csv",
    output_csv_path="output.csv"
)
```

## Advanced Usage

### Custom Action Approval Logic

```python
def custom_approval_logic(actions):
    """Custom approval function"""
    for action in actions:
        # Auto-approve high-confidence renamings
        if action.action_type == 'rename_column':
            if action.details['similarity_score'] > 0.9:
                action.approved = True

        # Always approve ontology mappings
        elif action.action_type == 'add_ontology_mappings':
            action.approved = True

        # Review value harmonizations manually
        elif action.action_type.startswith('harmonize'):
            # Custom review logic here
            pass

    return actions

# Use custom reviewer
reviewer = HumanReviewAgent(validator.data_dictionary)
approved = reviewer.review_actions(
    all_actions,
    interactive=False,
    custom_reviewer=custom_approval_logic
)
```

### Workflow History and Reporting

```python
# Process multiple files
for input_file in input_files:
    summary = orchestrator.process_file(
        input_csv_path=input_file,
        output_csv_path=input_file.replace('.csv', '_harmonized.csv')
    )

# Save complete history
orchestrator.save_workflow_history("workflow_history.json")

# Generate reports
for summary in orchestrator.workflow_history:
    report = orchestrator.generate_summary_report(summary)
    print(report)
```

### Agent Action History

```python
# Get action history from an agent
history = column_renamer.get_action_history()

# Save agent history
column_renamer.save_history("column_renamer_history.json")

# Load agent history
column_renamer.load_history("column_renamer_history.json")
```

## Examples

See the `examples/` directory for complete examples:

- `simple_example.py`: Basic validation and harmonization workflow
- `advanced_example.py`: Custom configurations, batch approval, learned mappings

Run examples:
```bash
# Run simple example
python examples/simple_example.py

# Run all advanced examples
python examples/advanced_example.py

# Run specific advanced example
python examples/advanced_example.py --example 3
```

## Output Files

### Harmonized CSV
- Contains corrected column names and values
- Includes original columns (suffixed with `_original`) for provenance
- Adds new ontology mapping columns (e.g., `disease_ontology`, `disease_ontology_label`)

### Validation Report
- Lists extra columns not in data dictionary
- Shows valid columns
- Details invalid columns with error messages
- Reports missing required columns

### Workflow History (JSON)
- Complete record of all actions proposed and applied
- Timestamps and duration
- Before/after validation statistics
- Full action details for reproducibility

## Architecture

```
metadata_validator/
├── __init__.py              # Package exports
├── validator.py             # Core validation engine
├── orchestrator.py          # Workflow coordinator
├── agents/
│   ├── __init__.py
│   ├── base_agent.py        # Base agent class
│   ├── column_renamer_agent.py
│   ├── value_harmonizer_agent.py
│   ├── ontology_matcher_agent.py
│   └── human_review_agent.py
├── examples/
│   ├── simple_example.py
│   └── advanced_example.py
├── requirements.txt
└── README.md
```

## Best Practices

1. **Always Keep Originals**: The framework preserves original values for provenance
2. **Review High-Impact Changes**: Use human review for value harmonizations
3. **Build Learned Mappings**: Export and reuse learned value mappings across datasets
4. **Validate Incrementally**: Run validation, review results, then run harmonization
5. **Track History**: Save workflow history for reproducibility and auditing

## Troubleshooting

### Import Errors

If you get import errors, make sure the parent directory is in your Python path:

```python
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
```

### Ontology Backend Issues

- **Pronto**: Install with `pip install pronto` and provide valid OBO file path
- **BioPortal**: Get API key from https://bioportal.bioontology.org/
- **Custom**: Use the default built-in mappings (no external dependencies)

### Performance Optimization

- Use auto-approve mode for large datasets in production pipelines
- Set higher similarity thresholds to reduce false positives
- Process files in batches
- Cache ontology lookups for repeated terms

## Contributing

To extend the framework:

1. Create new agents by inheriting from `BaseAgent`
2. Implement `analyze()` and `apply_actions()` methods
3. Register agent in `agents/__init__.py`
4. Add agent to orchestrator workflow

## License

This framework is part of the curatedMetagenomicDataCuration project.

## Citation

If you use this framework, please cite the curatedMetagenomicData project.
