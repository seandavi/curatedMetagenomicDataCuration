"""
Orchestrator for the Metadata Validation and Harmonization System

This module coordinates the execution of multiple agents to validate,
harmonize, and enrich metadata files.
"""

import pandas as pd
from pathlib import Path
from typing import Dict, List, Any, Optional
import json
from datetime import datetime

from .validator import MetadataValidator
from .agents import (
    ColumnRenamerAgent,
    ValueHarmonizerAgent,
    OntologyMatcherAgent
)
from .agents.human_review_agent import HumanReviewAgent
from .agents.base_agent import AgentAction


class MetadataOrchestrator:
    """
    Orchestrates the complete metadata validation and harmonization workflow.

    The workflow consists of:
    1. Validation - Check metadata against data dictionary
    2. Column Renaming - Suggest and apply column name corrections
    3. Value Harmonization - Fix incorrect values
    4. Ontology Matching - Map terms to ontology concepts
    5. Human Review - Allow human approval of changes (optional)
    """

    def __init__(
        self,
        data_dictionary_path: str,
        ontology_backend: str = 'custom',
        ontology_config: Optional[Dict[str, Any]] = None,
        human_review: bool = True,
        auto_approve: bool = False
    ):
        """
        Initialize the orchestrator.

        Args:
            data_dictionary_path: Path to the data dictionary CSV
            ontology_backend: Backend for ontology matching ('pronto', 'bioportal', 'custom')
            ontology_config: Configuration for ontology backend
            human_review: Whether to enable human review of actions
            auto_approve: If True, automatically approve all actions
        """
        # Initialize validator
        self.validator = MetadataValidator(data_dictionary_path)

        # Initialize agents
        self.column_renamer = ColumnRenamerAgent(
            self.validator.data_dictionary,
            similarity_threshold=0.6
        )

        self.value_harmonizer = ValueHarmonizerAgent(
            self.validator.data_dictionary,
            similarity_threshold=0.8
        )

        self.ontology_matcher = OntologyMatcherAgent(
            self.validator.data_dictionary,
            backend=ontology_backend,
            ontology_config=ontology_config
        )

        self.human_reviewer = HumanReviewAgent(
            self.validator.data_dictionary,
            auto_approve=auto_approve
        )

        self.human_review_enabled = human_review
        self.workflow_history: List[Dict[str, Any]] = []

    def process_file(
        self,
        input_csv_path: str,
        output_csv_path: Optional[str] = None,
        report_path: Optional[str] = None,
        interactive_review: bool = True
    ) -> Dict[str, Any]:
        """
        Process a metadata CSV file through the complete workflow.

        Args:
            input_csv_path: Path to input CSV file
            output_csv_path: Path to save harmonized output (optional)
            report_path: Path to save validation report (optional)
            interactive_review: Whether to use interactive human review

        Returns:
            Dictionary containing results and statistics
        """
        workflow_start = datetime.now()
        print(f"\n{'=' * 80}")
        print(f"METADATA HARMONIZATION WORKFLOW")
        print(f"Input file: {input_csv_path}")
        print(f"{'=' * 80}\n")

        # Step 1: Load data
        print("Step 1: Loading data...")
        df = pd.read_csv(input_csv_path)
        print(f"  Loaded {len(df)} rows, {len(df.columns)} columns\n")

        # Step 2: Initial validation
        print("Step 2: Validating against data dictionary...")
        validation_results = self.validator.validate_file(input_csv_path)
        print(f"  Extra columns: {len(validation_results['extra_columns'])}")
        print(f"  Valid columns: {len(validation_results['valid_columns'])}")
        print(f"  Invalid columns: {len(validation_results['invalid_columns'])}")
        print(f"  Missing required: {len(validation_results['missing_required_columns'])}\n")

        # Generate validation report
        if report_path:
            report = self.validator.generate_report(validation_results)
            with open(report_path, 'w') as f:
                f.write(report)
            print(f"  Validation report saved to: {report_path}\n")

        # Step 3: Collect proposed actions from all agents
        print("Step 3: Analyzing and proposing corrections...")
        all_actions = []

        # Column renaming suggestions
        print("  Running Column Renamer Agent...")
        rename_actions = self.column_renamer.analyze(df, validation_results)
        all_actions.extend(rename_actions)
        print(f"    Proposed {len(rename_actions)} column renamings")

        # Value harmonization suggestions
        print("  Running Value Harmonizer Agent...")
        harmonize_actions = self.value_harmonizer.analyze(df, validation_results)
        all_actions.extend(harmonize_actions)
        print(f"    Proposed {len(harmonize_actions)} value harmonizations")

        # Ontology matching suggestions
        print("  Running Ontology Matcher Agent...")
        ontology_actions = self.ontology_matcher.analyze(df, validation_results)
        all_actions.extend(ontology_actions)
        print(f"    Proposed {len(ontology_actions)} ontology mappings\n")

        print(f"  Total proposed actions: {len(all_actions)}\n")

        # Step 4: Human review (if enabled)
        if self.human_review_enabled and all_actions:
            print("Step 4: Human review of proposed actions...")
            all_actions = self.human_reviewer.review_actions(
                all_actions,
                interactive=interactive_review
            )
            approved_actions = [a for a in all_actions if a.approved]
            print(f"  {len(approved_actions)}/{len(all_actions)} actions approved\n")
        else:
            print("Step 4: Skipping human review (auto-approve enabled)...")
            for action in all_actions:
                action.approved = True
            approved_actions = all_actions
            print(f"  All {len(all_actions)} actions auto-approved\n")

        # Step 5: Apply approved actions
        print("Step 5: Applying approved actions...")
        df_harmonized = df.copy()

        # Apply in order: rename columns, harmonize values, add ontology mappings
        rename_approved = [a for a in approved_actions if a.action_type == 'rename_column']
        if rename_approved:
            print(f"  Applying {len(rename_approved)} column renamings...")
            df_harmonized = self.column_renamer.apply_actions(df_harmonized, rename_approved)

        harmonize_approved = [a for a in approved_actions if a.action_type.startswith('harmonize')]
        if harmonize_approved:
            print(f"  Applying {len(harmonize_approved)} value harmonizations...")
            df_harmonized = self.value_harmonizer.apply_actions(df_harmonized, harmonize_approved)

        ontology_approved = [a for a in approved_actions if a.action_type == 'add_ontology_mappings']
        if ontology_approved:
            print(f"  Applying {len(ontology_approved)} ontology mappings...")
            df_harmonized = self.ontology_matcher.apply_actions(df_harmonized, ontology_approved)

        print(f"  Result: {len(df_harmonized)} rows, {len(df_harmonized.columns)} columns\n")

        # Step 6: Re-validate
        print("Step 6: Re-validating harmonized data...")
        if output_csv_path:
            temp_output = output_csv_path.replace('.csv', '_temp.csv')
            df_harmonized.to_csv(temp_output, index=False)
            validation_results_final = self.validator.validate_file(temp_output)
            Path(temp_output).unlink()  # Remove temp file
        else:
            # Create temporary file for validation
            temp_output = '/tmp/metadata_temp.csv'
            df_harmonized.to_csv(temp_output, index=False)
            validation_results_final = self.validator.validate_file(temp_output)
            Path(temp_output).unlink()

        print(f"  Extra columns: {len(validation_results_final['extra_columns'])}")
        print(f"  Valid columns: {len(validation_results_final['valid_columns'])}")
        print(f"  Invalid columns: {len(validation_results_final['invalid_columns'])}")
        print(f"  Missing required: {len(validation_results_final['missing_required_columns'])}\n")

        # Step 7: Save output
        if output_csv_path:
            print(f"Step 7: Saving harmonized data to: {output_csv_path}")
            df_harmonized.to_csv(output_csv_path, index=False)
            print(f"  Saved successfully\n")

        # Create workflow summary
        workflow_end = datetime.now()
        workflow_duration = (workflow_end - workflow_start).total_seconds()

        summary = {
            'input_file': input_csv_path,
            'output_file': output_csv_path,
            'workflow_start': workflow_start.isoformat(),
            'workflow_end': workflow_end.isoformat(),
            'duration_seconds': workflow_duration,
            'initial_validation': {
                'extra_columns': len(validation_results['extra_columns']),
                'valid_columns': len(validation_results['valid_columns']),
                'invalid_columns': len(validation_results['invalid_columns']),
                'missing_required': len(validation_results['missing_required_columns'])
            },
            'final_validation': {
                'extra_columns': len(validation_results_final['extra_columns']),
                'valid_columns': len(validation_results_final['valid_columns']),
                'invalid_columns': len(validation_results_final['invalid_columns']),
                'missing_required': len(validation_results_final['missing_required_columns'])
            },
            'actions_proposed': len(all_actions),
            'actions_approved': len(approved_actions),
            'actions_applied': len([a for a in approved_actions if a.applied]),
            'all_actions': [a.to_dict() for a in all_actions]
        }

        self.workflow_history.append(summary)

        print(f"{'=' * 80}")
        print(f"WORKFLOW COMPLETE")
        print(f"Duration: {workflow_duration:.2f} seconds")
        print(f"{'=' * 80}\n")

        return summary

    def save_workflow_history(self, filepath: str):
        """Save complete workflow history to JSON file"""
        with open(filepath, 'w') as f:
            json.dump(self.workflow_history, f, indent=2)

    def generate_summary_report(self, workflow_summary: Dict[str, Any]) -> str:
        """Generate a human-readable summary report"""
        report = []
        report.append("=" * 80)
        report.append("WORKFLOW SUMMARY REPORT")
        report.append("=" * 80)
        report.append("")
        report.append(f"Input file: {workflow_summary['input_file']}")
        report.append(f"Output file: {workflow_summary.get('output_file', 'N/A')}")
        report.append(f"Duration: {workflow_summary['duration_seconds']:.2f} seconds")
        report.append("")

        report.append("INITIAL VALIDATION:")
        report.append("-" * 80)
        initial = workflow_summary['initial_validation']
        report.append(f"  Extra columns: {initial['extra_columns']}")
        report.append(f"  Valid columns: {initial['valid_columns']}")
        report.append(f"  Invalid columns: {initial['invalid_columns']}")
        report.append(f"  Missing required columns: {initial['missing_required']}")
        report.append("")

        report.append("ACTIONS:")
        report.append("-" * 80)
        report.append(f"  Proposed: {workflow_summary['actions_proposed']}")
        report.append(f"  Approved: {workflow_summary['actions_approved']}")
        report.append(f"  Applied: {workflow_summary['actions_applied']}")
        report.append("")

        report.append("FINAL VALIDATION:")
        report.append("-" * 80)
        final = workflow_summary['final_validation']
        report.append(f"  Extra columns: {final['extra_columns']}")
        report.append(f"  Valid columns: {final['valid_columns']}")
        report.append(f"  Invalid columns: {final['invalid_columns']}")
        report.append(f"  Missing required columns: {final['missing_required']}")
        report.append("")

        # Improvement summary
        extra_improved = initial['extra_columns'] - final['extra_columns']
        invalid_improved = initial['invalid_columns'] - final['invalid_columns']

        report.append("IMPROVEMENTS:")
        report.append("-" * 80)
        report.append(f"  Extra columns resolved: {extra_improved}")
        report.append(f"  Invalid columns fixed: {invalid_improved}")
        report.append("")

        report.append("=" * 80)
        return "\n".join(report)
