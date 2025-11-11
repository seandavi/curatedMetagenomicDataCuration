"""
Value Harmonizer Agent

This agent analyzes values in columns that don't match the data dictionary
and suggests corrections based on fuzzy matching and common transformation patterns.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Tuple
from difflib import get_close_matches, SequenceMatcher
import re
from .base_agent import BaseAgent, AgentAction


class ValueHarmonizerAgent(BaseAgent):
    """
    Agent that harmonizes incorrect values to match data dictionary specifications.

    Uses fuzzy matching, common transformation patterns, and learned mappings
    to suggest value corrections.
    """

    def __init__(self, data_dictionary: Dict[str, Any], similarity_threshold: float = 0.8):
        """
        Initialize the Value Harmonizer Agent.

        Args:
            data_dictionary: Dictionary mapping column names to DataDictionaryEntry objects
            similarity_threshold: Minimum similarity score (0-1) to suggest a match
        """
        super().__init__("ValueHarmonizerAgent", data_dictionary)
        self.similarity_threshold = similarity_threshold
        self.learned_mappings: Dict[str, Dict[str, str]] = {}  # Column -> {bad_value: good_value}

    def analyze(self, df: pd.DataFrame, validation_results: Optional[Dict] = None) -> List[AgentAction]:
        """
        Analyze invalid values and suggest corrections.

        Args:
            df: DataFrame to analyze
            validation_results: Results from MetadataValidator (should contain 'invalid_columns')

        Returns:
            List of proposed value harmonization actions
        """
        proposed_actions = []

        if validation_results and 'invalid_columns' in validation_results:
            invalid_columns = validation_results['invalid_columns']
        else:
            # Would need to validate if not provided
            return proposed_actions

        for col_result in invalid_columns:
            column_name = col_result.column_name

            if column_name not in self.data_dictionary:
                continue

            dict_entry = self.data_dictionary[column_name]

            # Get allowed values
            allowed_values = dict_entry.get_allowed_values_list()

            if not allowed_values:
                # For regex patterns, we can't suggest specific corrections easily
                # But we can identify common patterns
                action = self._analyze_regex_patterns(df, column_name, dict_entry, col_result)
                if action:
                    proposed_actions.append(action)
                continue

            # For enumerated values, suggest corrections
            value_mappings = self._suggest_value_mappings(
                df, column_name, allowed_values, col_result
            )

            if value_mappings:
                action = self.create_action(
                    action_type='harmonize_values',
                    description=f"Harmonize {len(value_mappings)} invalid values in '{column_name}'",
                    details={
                        'column': column_name,
                        'mappings': value_mappings,
                        'allowed_values': allowed_values,
                        'keep_original': True  # Keep original for provenance
                    }
                )
                proposed_actions.append(action)

        return proposed_actions

    def _suggest_value_mappings(
        self,
        df: pd.DataFrame,
        column_name: str,
        allowed_values: List[str],
        col_result: Any
    ) -> Dict[str, str]:
        """
        Suggest mappings from invalid values to valid values.

        Args:
            df: DataFrame containing the data
            column_name: Name of the column
            allowed_values: List of valid values from data dictionary
            col_result: ValidationResult for this column

        Returns:
            Dictionary mapping invalid values to suggested valid values
        """
        mappings = {}

        # Get unique invalid values
        invalid_values = set(col_result.invalid_values)

        for invalid_value in invalid_values:
            # Skip NA values
            if pd.isna(invalid_value) or invalid_value == '' or invalid_value == 'NA':
                continue

            # Check if we have a learned mapping
            if column_name in self.learned_mappings:
                if invalid_value in self.learned_mappings[column_name]:
                    mappings[invalid_value] = self.learned_mappings[column_name][invalid_value]
                    continue

            # Try to find a close match
            suggestion = self._find_best_match(invalid_value, allowed_values)
            if suggestion:
                mappings[invalid_value] = suggestion

        return mappings

    def _find_best_match(self, value: str, allowed_values: List[str]) -> Optional[str]:
        """
        Find the best matching allowed value for an invalid value.

        Args:
            value: Invalid value to match
            allowed_values: List of allowed values

        Returns:
            Best matching allowed value, or None if no good match
        """
        value_str = str(value)

        # Try exact match (case-insensitive)
        for allowed in allowed_values:
            if value_str.lower() == allowed.lower():
                return allowed

        # Try common transformations
        transformed = self._apply_common_transformations(value_str)
        for allowed in allowed_values:
            if transformed.lower() == allowed.lower():
                return allowed

        # Try fuzzy matching
        matches = get_close_matches(
            value_str,
            allowed_values,
            n=1,
            cutoff=self.similarity_threshold
        )
        if matches:
            return matches[0]

        # Try partial matching
        for allowed in allowed_values:
            similarity = SequenceMatcher(None, value_str.lower(), allowed.lower()).ratio()
            if similarity >= self.similarity_threshold:
                return allowed

        return None

    def _apply_common_transformations(self, value: str) -> str:
        """
        Apply common transformations that might fix the value.

        Examples:
        - "male" -> "Male"
        - "Type-2-Diabetes" -> "Type 2 Diabetes"
        - "healthy_control" -> "Healthy"
        """
        # Strip whitespace
        value = value.strip()

        # Replace underscores and hyphens with spaces
        value = re.sub(r'[-_]+', ' ', value)

        # Normalize multiple spaces
        value = re.sub(r'\s+', ' ', value)

        # Title case for common patterns
        if value.lower() in ['male', 'female', 'yes', 'no']:
            value = value.capitalize()

        return value

    def _analyze_regex_patterns(
        self,
        df: pd.DataFrame,
        column_name: str,
        dict_entry: Any,
        col_result: Any
    ) -> Optional[AgentAction]:
        """
        Analyze patterns in regex-validated columns.

        For regex patterns, we look for common issues like:
        - Extra whitespace
        - Wrong case
        - Missing/extra characters

        Args:
            df: DataFrame
            column_name: Column name
            dict_entry: Data dictionary entry
            col_result: Validation result

        Returns:
            Action if patterns are found, None otherwise
        """
        pattern = dict_entry.allowed_values
        invalid_values = set(col_result.invalid_values)

        transformations = {}
        for value in invalid_values:
            if pd.isna(value) or value == '':
                continue

            value_str = str(value).strip()

            # Try with whitespace removed
            no_space = re.sub(r'\s+', '', value_str)
            if re.match(f'^{pattern}$', no_space):
                transformations[value] = no_space
                continue

            # Try with case changes
            if re.match(f'^{pattern}$', value_str.upper()):
                transformations[value] = value_str.upper()
            elif re.match(f'^{pattern}$', value_str.lower()):
                transformations[value] = value_str.lower()

        if transformations:
            return self.create_action(
                action_type='harmonize_regex_values',
                description=f"Fix {len(transformations)} pattern violations in '{column_name}'",
                details={
                    'column': column_name,
                    'transformations': transformations,
                    'pattern': pattern,
                    'keep_original': True
                }
            )

        return None

    def apply_actions(self, df: pd.DataFrame, actions: List[AgentAction]) -> pd.DataFrame:
        """
        Apply approved value harmonization actions.

        Args:
            df: DataFrame to modify
            actions: List of approved harmonization actions

        Returns:
            Modified DataFrame with harmonized values
        """
        df_copy = df.copy()

        for action in actions:
            if not action.approved:
                continue

            if action.action_type == 'harmonize_values':
                column = action.details['column']
                mappings = action.details['mappings']
                keep_original = action.details.get('keep_original', True)

                if column in df_copy.columns:
                    if keep_original:
                        # Keep original values for provenance
                        df_copy[f"{column}_original"] = df_copy[column]

                    # Apply mappings
                    df_copy[column] = df_copy[column].replace(mappings)
                    action.applied = True

            elif action.action_type == 'harmonize_regex_values':
                column = action.details['column']
                transformations = action.details['transformations']
                keep_original = action.details.get('keep_original', True)

                if column in df_copy.columns:
                    if keep_original:
                        df_copy[f"{column}_original"] = df_copy[column]

                    df_copy[column] = df_copy[column].replace(transformations)
                    action.applied = True

        return df_copy

    def learn_mapping(self, column: str, incorrect_value: str, correct_value: str):
        """
        Teach the agent a new value mapping.

        Args:
            column: Column name
            incorrect_value: The incorrect value
            correct_value: The correct value to map to
        """
        if column not in self.learned_mappings:
            self.learned_mappings[column] = {}

        self.learned_mappings[column][incorrect_value] = correct_value

    def export_learned_mappings(self, filepath: str):
        """Export learned mappings to a JSON file"""
        import json
        with open(filepath, 'w') as f:
            json.dump(self.learned_mappings, f, indent=2)

    def import_learned_mappings(self, filepath: str):
        """Import learned mappings from a JSON file"""
        import json
        with open(filepath, 'r') as f:
            self.learned_mappings = json.load(f)
