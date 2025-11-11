#!/usr/bin/env python3
"""
Advanced Example: Custom Configuration and Agent Usage

This example demonstrates advanced features:
- Custom ontology mappings
- Manual agent usage
- Learned value mappings
- Batch approval strategies
"""

import sys
from pathlib import Path
import json

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from metadata_validator import (
    MetadataValidator,
    ColumnRenamerAgent,
    ValueHarmonizerAgent,
    OntologyMatcherAgent,
    HumanReviewAgent
)
import pandas as pd


def example_1_manual_validation():
    """Example 1: Manual validation without orchestrator"""
    print("=" * 80)
    print("EXAMPLE 1: Manual Validation")
    print("=" * 80)

    # Initialize validator
    validator = MetadataValidator("../inst/extdata/cMD_data_dictionary.csv")

    # Validate a file
    results = validator.validate_file("sample_metadata.csv")

    # Generate and print report
    report = validator.generate_report(results)
    print(report)


def example_2_custom_agents():
    """Example 2: Using agents independently with custom configuration"""
    print("\n" + "=" * 80)
    print("EXAMPLE 2: Custom Agent Usage")
    print("=" * 80)

    # Load data dictionary
    validator = MetadataValidator("../inst/extdata/cMD_data_dictionary.csv")
    data_dict = validator.data_dictionary

    # Load data
    df = pd.read_csv("sample_metadata.csv")

    # Validate
    validation_results = validator.validate_file("sample_metadata.csv")

    # Initialize agents with custom settings
    column_renamer = ColumnRenamerAgent(
        data_dict,
        similarity_threshold=0.5  # Lower threshold for more suggestions
    )

    value_harmonizer = ValueHarmonizerAgent(
        data_dict,
        similarity_threshold=0.7
    )

    # Custom ontology configuration
    custom_ontology_mappings = {
        'Custom Disease': {
            'ontology': 'CUSTOM',
            'term_id': 'CUSTOM:001',
            'term_label': 'Custom Disease Name',
            'definition': 'A custom disease definition'
        }
    }

    ontology_matcher = OntologyMatcherAgent(
        data_dict,
        backend='custom',
        ontology_config={'mappings': custom_ontology_mappings}
    )

    # Analyze with each agent
    print("\nColumn Renamer Agent:")
    rename_actions = column_renamer.analyze(df, validation_results)
    for action in rename_actions:
        print(f"  {action.description}")

    print("\nValue Harmonizer Agent:")
    harmonize_actions = value_harmonizer.analyze(df, validation_results)
    for action in harmonize_actions:
        print(f"  {action.description}")

    print("\nOntology Matcher Agent:")
    ontology_actions = ontology_matcher.analyze(df, validation_results)
    for action in ontology_actions:
        print(f"  {action.description}")


def example_3_learned_mappings():
    """Example 3: Teaching the harmonizer custom mappings"""
    print("\n" + "=" * 80)
    print("EXAMPLE 3: Learned Value Mappings")
    print("=" * 80)

    validator = MetadataValidator("../inst/extdata/cMD_data_dictionary.csv")
    harmonizer = ValueHarmonizerAgent(validator.data_dictionary)

    # Teach custom mappings
    harmonizer.learn_mapping('sex', 'M', 'Male')
    harmonizer.learn_mapping('sex', 'F', 'Female')
    harmonizer.learn_mapping('sex', 'm', 'Male')
    harmonizer.learn_mapping('sex', 'f', 'Female')
    harmonizer.learn_mapping('disease', 'T2D', 'Type 2 Diabetes Mellitus')
    harmonizer.learn_mapping('disease', 'T1D', 'Type 1 Diabetes Mellitus')
    harmonizer.learn_mapping('disease', 'CRC', 'Colorectal Carcinoma')

    print("Learned mappings:")
    print(json.dumps(harmonizer.learned_mappings, indent=2))

    # Export for reuse
    harmonizer.export_learned_mappings("learned_mappings.json")
    print("\nLearned mappings exported to: learned_mappings.json")


def example_4_batch_approval():
    """Example 4: Batch approval strategies"""
    print("\n" + "=" * 80)
    print("EXAMPLE 4: Batch Approval Strategies")
    print("=" * 80)

    validator = MetadataValidator("../inst/extdata/cMD_data_dictionary.csv")
    df = pd.read_csv("sample_metadata.csv")
    validation_results = validator.validate_file("sample_metadata.csv")

    # Create agents
    renamer = ColumnRenamerAgent(validator.data_dictionary)
    harmonizer = ValueHarmonizerAgent(validator.data_dictionary)
    ontology = OntologyMatcherAgent(validator.data_dictionary)

    # Collect all actions
    all_actions = []
    all_actions.extend(renamer.analyze(df, validation_results))
    all_actions.extend(harmonizer.analyze(df, validation_results))
    all_actions.extend(ontology.analyze(df, validation_results))

    # Create review agent
    reviewer = HumanReviewAgent(validator.data_dictionary)

    # Strategy 1: Approve all ontology mappings (safe operation)
    print("\nStrategy 1: Auto-approve all ontology mappings")
    reviewer.batch_approve_by_type(all_actions, 'add_ontology_mappings')

    # Strategy 2: Approve all actions from a trusted agent
    print("Strategy 2: Auto-approve all column renaming with high confidence")
    for action in all_actions:
        if action.action_type == 'rename_column':
            if action.details['similarity_score'] > 0.9:
                action.approved = True

    # Strategy 3: Review only critical changes interactively
    print("Strategy 3: Interactive review for value harmonization only")
    critical_actions = [a for a in all_actions if a.action_type.startswith('harmonize')]
    if critical_actions:
        reviewer.review_actions(critical_actions, interactive=True)

    # Generate approval report
    approval_report = reviewer.generate_approval_report(all_actions)
    print("\n" + approval_report)


def example_5_ontology_backends():
    """Example 5: Using different ontology backends"""
    print("\n" + "=" * 80)
    print("EXAMPLE 5: Ontology Backend Options")
    print("=" * 80)

    validator = MetadataValidator("../inst/extdata/cMD_data_dictionary.csv")

    # Option 1: Custom backend (default, no external dependencies)
    print("\nOption 1: Custom backend (built-in mappings)")
    matcher_custom = OntologyMatcherAgent(
        validator.data_dictionary,
        backend='custom'
    )
    print("  Status: Ready to use")

    # Option 2: Pronto backend (requires pronto library and OBO file)
    print("\nOption 2: Pronto backend (local OBO ontology)")
    print("  Requirements: pip install pronto")
    print("  Usage:")
    print("    matcher_pronto = OntologyMatcherAgent(")
    print("        data_dictionary,")
    print("        backend='pronto',")
    print("        ontology_config={'ontology_path': '/path/to/ontology.obo'}")
    print("    )")

    # Option 3: BioPortal backend (requires API key)
    print("\nOption 3: BioPortal API backend")
    print("  Requirements: pip install requests")
    print("  Get API key from: https://bioportal.bioontology.org/")
    print("  Usage:")
    print("    matcher_bioportal = OntologyMatcherAgent(")
    print("        data_dictionary,")
    print("        backend='bioportal',")
    print("        ontology_config={")
    print("            'api_key': 'your-api-key',")
    print("            'ontologies': ['NCIT', 'DOID', 'HP']")
    print("        }")
    print("    )")


def main():
    """Run all examples"""
    import argparse

    parser = argparse.ArgumentParser(description="Advanced metadata validator examples")
    parser.add_argument(
        '--example',
        type=int,
        choices=[1, 2, 3, 4, 5],
        help="Run a specific example (1-5), or run all if not specified"
    )

    args = parser.parse_args()

    examples = {
        1: example_1_manual_validation,
        2: example_2_custom_agents,
        3: example_3_learned_mappings,
        4: example_4_batch_approval,
        5: example_5_ontology_backends
    }

    if args.example:
        examples[args.example]()
    else:
        for example_func in examples.values():
            try:
                example_func()
            except FileNotFoundError:
                print(f"\nSkipping example (sample file not found)")
            except Exception as e:
                print(f"\nError in example: {e}")


if __name__ == "__main__":
    main()
