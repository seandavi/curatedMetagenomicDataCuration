"""
Core Validator Module for Curated Metagenomic Data Curation

This module provides functionality to validate CSV metadata files against
a data dictionary with support for regex patterns and enumerated values.
"""

import pandas as pd
import re
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class ValidationResult:
    """Results from validating a single column"""
    column_name: str
    is_valid: bool
    invalid_rows: List[int] = field(default_factory=list)
    invalid_values: List[str] = field(default_factory=list)
    error_messages: List[str] = field(default_factory=list)
    allowed_values: Optional[str] = None
    column_class: Optional[str] = None


@dataclass
class DataDictionaryEntry:
    """Represents a single entry in the data dictionary"""
    col_name: str
    col_class: str
    unique: str
    required: str
    multiple_values: bool
    description: str
    allowed_values: str
    delimiter: str
    separator: str
    dynamic_enum: str
    dynamic_enum_property: str

    def is_regex_pattern(self) -> bool:
        """Check if allowed_values contains a regex pattern"""
        if pd.isna(self.allowed_values) or not self.allowed_values or self.allowed_values == 'NA':
            return False
        # Check if it's not a simple pipe-separated list
        # Regex patterns typically contain special characters like [], ^, $, +, etc.
        regex_chars = set('[]()+*.^$\\')
        return any(char in str(self.allowed_values) for char in regex_chars)

    def get_allowed_values_list(self) -> Optional[List[str]]:
        """Get list of allowed values if it's an enumeration"""
        if pd.isna(self.allowed_values) or not self.allowed_values or self.allowed_values == 'NA' or self.is_regex_pattern():
            return None
        return [v.strip() for v in str(self.allowed_values).split('|')]

    def validate_value(self, value: str) -> bool:
        """Validate a single value against this dictionary entry"""
        if pd.isna(value) or value == '' or value == 'NA':
            return self.required != 'required'

        value_str = str(value)

        # Handle multiple values
        if self.multiple_values and self.delimiter != 'NA':
            values = [v.strip() for v in value_str.split(self.delimiter)]
            return all(self._validate_single_value(v) for v in values)

        return self._validate_single_value(value_str)

    def _validate_single_value(self, value: str) -> bool:
        """Validate a single value (not split by delimiter)"""
        if pd.isna(self.allowed_values) or not self.allowed_values or self.allowed_values == 'NA':
            return True

        if self.is_regex_pattern():
            # Validate against regex
            try:
                pattern = re.compile(f'^{str(self.allowed_values)}$')
                return bool(pattern.match(value))
            except re.error:
                # If regex is invalid, treat as exact match
                return value == str(self.allowed_values)
        else:
            # Validate against enumerated values
            allowed = self.get_allowed_values_list()
            return allowed is None or value in allowed


class MetadataValidator:
    """
    Validates metadata CSV files against a data dictionary.

    The validator checks:
    - Extra columns not in the data dictionary
    - Columns that validate correctly
    - Columns with invalid values (with details about what failed)
    """

    def __init__(self, data_dictionary_path: str):
        """
        Initialize validator with a data dictionary.

        Args:
            data_dictionary_path: Path to the data dictionary CSV file
        """
        self.data_dictionary_path = Path(data_dictionary_path)
        self.data_dictionary: Dict[str, DataDictionaryEntry] = {}
        self._load_data_dictionary()

    def _load_data_dictionary(self):
        """Load and parse the data dictionary"""
        df = pd.read_csv(self.data_dictionary_path)

        for _, row in df.iterrows():
            entry = DataDictionaryEntry(
                col_name=row['ColName'],
                col_class=row['ColClass'],
                unique=row['Unique'],
                required=row['Required'],
                multiple_values=row['MultipleValues'],
                description=row['Description'],
                allowed_values=row['AllowedValues'],
                delimiter=row['Delimiter'],
                separator=row['Separater'],
                dynamic_enum=row['DynamicEnum'],
                dynamic_enum_property=row['DynamicEnumProperty']
            )
            self.data_dictionary[entry.col_name] = entry

    def validate_file(self, csv_path: str) -> Dict[str, any]:
        """
        Validate a CSV file against the data dictionary.

        Args:
            csv_path: Path to the CSV file to validate

        Returns:
            Dictionary containing validation results with keys:
            - 'extra_columns': Columns not in data dictionary
            - 'valid_columns': Columns that validate successfully
            - 'invalid_columns': Columns with validation errors
            - 'missing_required_columns': Required columns that are missing
        """
        df = pd.read_csv(csv_path)

        results = {
            'extra_columns': [],
            'valid_columns': [],
            'invalid_columns': [],
            'missing_required_columns': []
        }

        # Check for extra columns
        for col in df.columns:
            if col not in self.data_dictionary:
                results['extra_columns'].append(col)

        # Check for missing required columns
        for col_name, entry in self.data_dictionary.items():
            if entry.required == 'required' and col_name not in df.columns:
                results['missing_required_columns'].append(col_name)

        # Validate columns that exist in both
        for col in df.columns:
            if col in self.data_dictionary:
                validation_result = self._validate_column(df, col)
                if validation_result.is_valid:
                    results['valid_columns'].append(validation_result)
                else:
                    results['invalid_columns'].append(validation_result)

        return results

    def _validate_column(self, df: pd.DataFrame, column_name: str) -> ValidationResult:
        """
        Validate a single column in the dataframe.

        Args:
            df: The dataframe containing the data
            column_name: Name of the column to validate

        Returns:
            ValidationResult object with validation details
        """
        entry = self.data_dictionary[column_name]
        result = ValidationResult(
            column_name=column_name,
            is_valid=True,
            allowed_values=entry.allowed_values,
            column_class=entry.col_class
        )

        for idx, value in enumerate(df[column_name]):
            if not entry.validate_value(value):
                result.is_valid = False
                result.invalid_rows.append(idx)
                result.invalid_values.append(str(value))

                # Create helpful error message
                if entry.is_regex_pattern():
                    error_msg = f"Row {idx}: '{value}' does not match pattern '{entry.allowed_values}'"
                else:
                    allowed = entry.get_allowed_values_list()
                    if allowed:
                        error_msg = f"Row {idx}: '{value}' not in allowed values: {allowed[:5]}..."
                    else:
                        error_msg = f"Row {idx}: '{value}' is invalid"

                result.error_messages.append(error_msg)

        return result

    def generate_report(self, validation_results: Dict[str, any]) -> str:
        """
        Generate a human-readable validation report.

        Args:
            validation_results: Results from validate_file()

        Returns:
            Formatted string report
        """
        report = []
        report.append("=" * 80)
        report.append("METADATA VALIDATION REPORT")
        report.append("=" * 80)
        report.append("")

        # Extra columns
        if validation_results['extra_columns']:
            report.append("EXTRA COLUMNS (not in data dictionary):")
            report.append("-" * 80)
            for col in validation_results['extra_columns']:
                report.append(f"  - {col}")
            report.append("")
        else:
            report.append("✓ No extra columns found")
            report.append("")

        # Missing required columns
        if validation_results['missing_required_columns']:
            report.append("MISSING REQUIRED COLUMNS:")
            report.append("-" * 80)
            for col in validation_results['missing_required_columns']:
                report.append(f"  - {col}")
            report.append("")
        else:
            report.append("✓ All required columns present")
            report.append("")

        # Valid columns
        report.append(f"VALID COLUMNS ({len(validation_results['valid_columns'])}):")
        report.append("-" * 80)
        for result in validation_results['valid_columns']:
            report.append(f"  ✓ {result.column_name}")
        report.append("")

        # Invalid columns
        if validation_results['invalid_columns']:
            report.append(f"INVALID COLUMNS ({len(validation_results['invalid_columns'])}):")
            report.append("-" * 80)
            for result in validation_results['invalid_columns']:
                report.append(f"  ✗ {result.column_name}")
                report.append(f"    Invalid rows: {len(result.invalid_rows)}")
                report.append(f"    Allowed values: {result.allowed_values[:100]}...")
                report.append(f"    Sample errors (showing first 5):")
                for error in result.error_messages[:5]:
                    report.append(f"      {error}")
                if len(result.error_messages) > 5:
                    report.append(f"      ... and {len(result.error_messages) - 5} more errors")
                report.append("")
        else:
            report.append("✓ All columns validated successfully!")
            report.append("")

        report.append("=" * 80)
        return "\n".join(report)
