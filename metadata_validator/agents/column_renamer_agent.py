"""
Column Renamer Agent

This agent analyzes column names that don't match the data dictionary
and suggests corrections based on similarity matching.
"""

import pandas as pd
from typing import Dict, List, Any, Optional, Tuple
from difflib import SequenceMatcher
import re
from .base_agent import BaseAgent, AgentAction


class ColumnRenamerAgent(BaseAgent):
    """
    Agent that suggests column name corrections to match the data dictionary.

    Uses fuzzy string matching to find the most likely correct column name
    for columns that don't exist in the data dictionary.
    """

    def __init__(self, data_dictionary: Dict[str, Any], similarity_threshold: float = 0.6):
        """
        Initialize the Column Renamer Agent.

        Args:
            data_dictionary: Dictionary mapping column names to DataDictionaryEntry objects
            similarity_threshold: Minimum similarity score (0-1) to suggest a match
        """
        super().__init__("ColumnRenamerAgent", data_dictionary)
        self.similarity_threshold = similarity_threshold

    def analyze(self, df: pd.DataFrame, validation_results: Optional[Dict] = None) -> List[AgentAction]:
        """
        Analyze extra columns and suggest renamings.

        Args:
            df: DataFrame to analyze
            validation_results: Results from MetadataValidator (should contain 'extra_columns')

        Returns:
            List of proposed renaming actions
        """
        proposed_actions = []

        if validation_results and 'extra_columns' in validation_results:
            extra_columns = validation_results['extra_columns']
        else:
            # Find extra columns if not provided
            extra_columns = [col for col in df.columns if col not in self.data_dictionary]

        for column in extra_columns:
            suggestions = self._find_similar_columns(column)

            if suggestions:
                best_match, score = suggestions[0]

                action = self.create_action(
                    action_type='rename_column',
                    description=f"Rename '{column}' to '{best_match}' (similarity: {score:.2f})",
                    details={
                        'original_column': column,
                        'suggested_column': best_match,
                        'similarity_score': score,
                        'all_suggestions': suggestions,
                        'keep_original': True  # Keep original for provenance
                    }
                )
                proposed_actions.append(action)

        return proposed_actions

    def _find_similar_columns(self, column_name: str) -> List[Tuple[str, float]]:
        """
        Find similar column names in the data dictionary.

        Args:
            column_name: Name of the column to match

        Returns:
            List of (column_name, similarity_score) tuples, sorted by score descending
        """
        similarities = []

        for dict_column in self.data_dictionary.keys():
            score = self._calculate_similarity(column_name, dict_column)
            if score >= self.similarity_threshold:
                similarities.append((dict_column, score))

        # Sort by similarity score (highest first)
        similarities.sort(key=lambda x: x[1], reverse=True)
        return similarities

    def _calculate_similarity(self, str1: str, str2: str) -> float:
        """
        Calculate similarity between two strings using multiple methods.

        Uses a combination of:
        - Sequence matching
        - Case-insensitive comparison
        - Special handling for common separators (_, -, space)

        Args:
            str1: First string
            str2: Second string

        Returns:
            Similarity score between 0 and 1
        """
        # Normalize strings: lowercase and replace separators
        norm1 = self._normalize_string(str1)
        norm2 = self._normalize_string(str2)

        # Calculate sequence match ratio
        seq_ratio = SequenceMatcher(None, norm1, norm2).ratio()

        # Bonus for exact match after normalization
        if norm1 == norm2:
            return 1.0

        # Bonus for substring matches
        substring_bonus = 0
        if norm1 in norm2 or norm2 in norm1:
            substring_bonus = 0.2

        return min(1.0, seq_ratio + substring_bonus)

    def _normalize_string(self, s: str) -> str:
        """Normalize string for comparison"""
        # Convert to lowercase
        s = s.lower()
        # Replace common separators with a single character
        s = re.sub(r'[-_\s]+', '_', s)
        # Remove non-alphanumeric characters except underscore
        s = re.sub(r'[^a-z0-9_]', '', s)
        return s

    def apply_actions(self, df: pd.DataFrame, actions: List[AgentAction]) -> pd.DataFrame:
        """
        Apply approved column renaming actions.

        Args:
            df: DataFrame to modify
            actions: List of approved renaming actions

        Returns:
            Modified DataFrame with renamed columns
        """
        df_copy = df.copy()

        for action in actions:
            if action.action_type == 'rename_column' and action.approved:
                original = action.details['original_column']
                new_name = action.details['suggested_column']
                keep_original = action.details.get('keep_original', True)

                if original in df_copy.columns:
                    if keep_original:
                        # Keep original column for provenance
                        df_copy[f"{original}_original"] = df_copy[original]

                    # Rename the column
                    df_copy = df_copy.rename(columns={original: new_name})
                    action.applied = True

        return df_copy

    def suggest_manual_mapping(self, extra_columns: List[str]) -> Dict[str, List[str]]:
        """
        Generate a manual mapping guide for columns with no good matches.

        Args:
            extra_columns: List of column names not in data dictionary

        Returns:
            Dictionary mapping original columns to top suggestions
        """
        mapping_guide = {}

        for column in extra_columns:
            suggestions = self._find_similar_columns(column)
            if suggestions:
                mapping_guide[column] = [name for name, score in suggestions[:5]]
            else:
                mapping_guide[column] = ["No suggestions - manual review needed"]

        return mapping_guide
