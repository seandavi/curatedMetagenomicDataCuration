#!/usr/bin/env python3
"""
Simple Example: Basic Metadata Validation and Harmonization

This example demonstrates the basic workflow for validating and harmonizing
metadata files using the MetadataOrchestrator.
"""

import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from metadata_validator import MetadataOrchestrator


def main():
    # Path to data dictionary
    data_dict_path = "../inst/extdata/cMD_data_dictionary.csv"

    # Create orchestrator
    # - human_review=True: Enable human review of proposed actions
    # - auto_approve=False: Require manual approval for each action
    orchestrator = MetadataOrchestrator(
        data_dictionary_path=data_dict_path,
        ontology_backend='custom',  # Use built-in ontology mappings
        human_review=True,
        auto_approve=False
    )

    # Process a metadata file
    input_file = "sample_metadata.csv"  # Replace with your file
    output_file = "sample_metadata_harmonized.csv"
    report_file = "validation_report.txt"

    print(f"Processing metadata file: {input_file}")
    print(f"Output will be saved to: {output_file}")
    print()

    # Run the complete workflow
    summary = orchestrator.process_file(
        input_csv_path=input_file,
        output_csv_path=output_file,
        report_path=report_file,
        interactive_review=True  # Enable interactive terminal review
    )

    # Generate and print summary
    summary_report = orchestrator.generate_summary_report(summary)
    print(summary_report)

    # Save workflow history
    orchestrator.save_workflow_history("workflow_history.json")
    print("\nWorkflow history saved to: workflow_history.json")


if __name__ == "__main__":
    main()
