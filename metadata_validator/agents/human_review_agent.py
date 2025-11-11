"""
Human-in-the-Loop Review Agent

This agent facilitates human review and approval of proposed actions
from other agents before they are applied.
"""

import pandas as pd
from typing import Dict, List, Any, Optional, Callable
from .base_agent import BaseAgent, AgentAction


class HumanReviewAgent(BaseAgent):
    """
    Agent that facilitates human review of proposed actions.

    Provides interactive review interface for approving/rejecting
    actions proposed by other agents.
    """

    def __init__(self, data_dictionary: Dict[str, Any], auto_approve: bool = False):
        """
        Initialize the Human Review Agent.

        Args:
            data_dictionary: Dictionary mapping column names to DataDictionaryEntry objects
            auto_approve: If True, automatically approve all actions (for testing)
        """
        super().__init__("HumanReviewAgent", data_dictionary)
        self.auto_approve = auto_approve
        self.reviewed_actions: List[AgentAction] = []

    def analyze(self, df: pd.DataFrame, validation_results: Optional[Dict] = None) -> List[AgentAction]:
        """
        This agent doesn't generate new actions, only reviews existing ones.

        Args:
            df: DataFrame (not used)
            validation_results: Validation results (not used)

        Returns:
            Empty list
        """
        return []

    def review_actions(
        self,
        actions: List[AgentAction],
        interactive: bool = True,
        custom_reviewer: Optional[Callable] = None
    ) -> List[AgentAction]:
        """
        Review and approve/reject actions.

        Args:
            actions: List of actions to review
            interactive: If True, prompt user for approval
            custom_reviewer: Optional custom review function

        Returns:
            List of approved actions
        """
        if self.auto_approve:
            for action in actions:
                action.approved = True
            return actions

        if custom_reviewer:
            return custom_reviewer(actions)

        if interactive:
            return self._interactive_review(actions)

        return actions

    def _interactive_review(self, actions: List[AgentAction]) -> List[AgentAction]:
        """
        Interactive terminal-based review of actions.

        Args:
            actions: List of actions to review

        Returns:
            List of actions (with approval status set)
        """
        print("\n" + "=" * 80)
        print("HUMAN REVIEW - Please review the following proposed actions")
        print("=" * 80)

        approved_count = 0

        for idx, action in enumerate(actions, 1):
            print(f"\nAction {idx}/{len(actions)}")
            print("-" * 80)
            print(f"Agent: {action.agent_name}")
            print(f"Type: {action.action_type}")
            print(f"Description: {action.description}")
            print(f"\nDetails:")
            self._print_action_details(action)

            while True:
                response = input("\nApprove this action? (y/n/s/q): ").lower().strip()

                if response == 'y':
                    action.approved = True
                    approved_count += 1
                    print("✓ Approved")
                    break
                elif response == 'n':
                    action.approved = False
                    print("✗ Rejected")
                    break
                elif response == 's':
                    action.approved = False
                    print("⊘ Skipped")
                    break
                elif response == 'q':
                    print("\nReview cancelled. No more actions will be reviewed.")
                    return actions
                else:
                    print("Invalid input. Please enter 'y' (yes), 'n' (no), 's' (skip), or 'q' (quit)")

        print("\n" + "=" * 80)
        print(f"Review complete: {approved_count}/{len(actions)} actions approved")
        print("=" * 80 + "\n")

        self.reviewed_actions.extend(actions)
        return actions

    def _print_action_details(self, action: AgentAction):
        """Print action details in a formatted way"""
        details = action.details

        if action.action_type == 'rename_column':
            print(f"  Original column: '{details['original_column']}'")
            print(f"  Suggested column: '{details['suggested_column']}'")
            print(f"  Similarity score: {details['similarity_score']:.2f}")
            if details.get('all_suggestions'):
                print(f"  Other suggestions:")
                for name, score in details['all_suggestions'][1:4]:  # Show top 3 alternatives
                    print(f"    - {name} (score: {score:.2f})")

        elif action.action_type == 'harmonize_values':
            print(f"  Column: '{details['column']}'")
            print(f"  Number of mappings: {len(details['mappings'])}")
            print(f"  Sample mappings (showing first 5):")
            for i, (old_val, new_val) in enumerate(list(details['mappings'].items())[:5], 1):
                print(f"    {i}. '{old_val}' -> '{new_val}'")
            if len(details['mappings']) > 5:
                print(f"    ... and {len(details['mappings']) - 5} more")

        elif action.action_type == 'harmonize_regex_values':
            print(f"  Column: '{details['column']}'")
            print(f"  Pattern: {details['pattern']}")
            print(f"  Number of transformations: {len(details['transformations'])}")
            print(f"  Sample transformations (showing first 5):")
            for i, (old_val, new_val) in enumerate(list(details['transformations'].items())[:5], 1):
                print(f"    {i}. '{old_val}' -> '{new_val}'")

        elif action.action_type == 'add_ontology_mappings':
            print(f"  Column: '{details['column']}'")
            print(f"  New column: '{details['new_column_name']}'")
            print(f"  Number of matches: {len(details['matches'])}")
            print(f"  Sample matches (showing first 5):")
            for i, (term, matches) in enumerate(list(details['matches'].items())[:5], 1):
                match = matches[0] if matches else {}
                term_id = match.get('term_id', 'N/A')
                term_label = match.get('term_label', 'N/A')
                ontology = match.get('ontology', 'N/A')
                print(f"    {i}. '{term}' -> {term_label} ({term_id}) [{ontology}]")

        else:
            # Generic detail printing
            for key, value in details.items():
                if isinstance(value, dict) and len(value) > 5:
                    print(f"  {key}: {len(value)} items")
                elif isinstance(value, list) and len(value) > 5:
                    print(f"  {key}: {len(value)} items")
                else:
                    print(f"  {key}: {value}")

    def apply_actions(self, df: pd.DataFrame, actions: List[AgentAction]) -> pd.DataFrame:
        """
        This agent doesn't apply actions directly.

        Args:
            df: DataFrame
            actions: Actions (not used)

        Returns:
            Unchanged DataFrame
        """
        return df

    def generate_approval_report(self, actions: List[AgentAction]) -> str:
        """
        Generate a report of approval decisions.

        Args:
            actions: List of reviewed actions

        Returns:
            Formatted report string
        """
        report = []
        report.append("=" * 80)
        report.append("HUMAN REVIEW REPORT")
        report.append("=" * 80)
        report.append("")

        approved = [a for a in actions if a.approved]
        rejected = [a for a in actions if not a.approved]

        report.append(f"Total actions reviewed: {len(actions)}")
        report.append(f"Approved: {len(approved)}")
        report.append(f"Rejected: {len(rejected)}")
        report.append("")

        if approved:
            report.append("APPROVED ACTIONS:")
            report.append("-" * 80)
            for action in approved:
                report.append(f"  ✓ [{action.agent_name}] {action.description}")
            report.append("")

        if rejected:
            report.append("REJECTED ACTIONS:")
            report.append("-" * 80)
            for action in rejected:
                report.append(f"  ✗ [{action.agent_name}] {action.description}")
            report.append("")

        report.append("=" * 80)
        return "\n".join(report)

    def batch_approve_by_agent(self, actions: List[AgentAction], agent_name: str):
        """
        Approve all actions from a specific agent.

        Args:
            actions: List of actions
            agent_name: Name of agent whose actions to approve
        """
        for action in actions:
            if action.agent_name == agent_name:
                action.approved = True

    def batch_approve_by_type(self, actions: List[AgentAction], action_type: str):
        """
        Approve all actions of a specific type.

        Args:
            actions: List of actions
            action_type: Type of actions to approve
        """
        for action in actions:
            if action.action_type == action_type:
                action.approved = True

    def batch_reject_by_agent(self, actions: List[AgentAction], agent_name: str):
        """
        Reject all actions from a specific agent.

        Args:
            actions: List of actions
            agent_name: Name of agent whose actions to reject
        """
        for action in actions:
            if action.agent_name == agent_name:
                action.approved = False
