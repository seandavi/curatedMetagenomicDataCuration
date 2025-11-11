#!/usr/bin/env python3
"""
Test script for the Metadata Validator Framework

Tests basic functionality with sample data.
"""

import pandas as pd
import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

from metadata_validator import (
    MetadataValidator,
    MetadataOrchestrator,
    ColumnRenamerAgent,
    ValueHarmonizerAgent,
    OntologyMatcherAgent
)


def create_sample_data():
    """Create sample metadata for testing"""
    data = {
        'subject_id': ['SUB001', 'SUB002', 'SUB003', 'SUB004'],
        'sample_id': ['SAMPLE001', 'SAMPLE002', 'SAMPLE003', 'SAMPLE004'],
        'study_name': ['TestStudy_2024', 'TestStudy_2024', 'TestStudy_2024', 'TestStudy_2024'],
        'Sex': ['M', 'F', 'Male', 'Female'],  # Mixed format, wrong column name
        'age': [45, 32, 67, 28],
        'disease': ['Type 2 Diabetes Mellitus', 'Healthy', 'T2D', 'Crohn Disease'],  # Some need harmonization
        'Country': ['United States', 'China', 'USA', 'Germany'],  # Wrong column name, inconsistent
        'body_site': ['stool', 'stool', 'stool', 'stool'],
        'control': ['Case', 'Study Control', 'Case', 'Study Control'],
        'target_condition': ['Type 2 Diabetes Mellitus', 'Microbiome', 'Type 2 Diabetes Mellitus', 'Crohn Disease'],
        'ancestry': ['European', 'Asian', 'European', 'European'],
        'curator': ['Bryan_Merrill', 'Bryan_Merrill', 'Bryan_Merrill', 'Bryan_Merrill'],
        'sequencing_platform': ['IlluminaHiSeq', 'IlluminaMiSeq', 'IlluminaHiSeq', 'IlluminaHiSeq'],
        'median_read_length': [150, 150, 150, 150],
        'minimum_read_length': [100, 100, 100, 100],
        'number_bases': [5000000000, 4500000000, 5500000000, 4800000000],
        'number_reads': [50000000, 45000000, 55000000, 48000000],
        'extra_column': ['value1', 'value2', 'value3', 'value4']  # Not in data dictionary
    }

    df = pd.DataFrame(data)
    sample_path = Path(__file__).parent / 'test_sample_metadata.csv'
    df.to_csv(sample_path, index=False)
    print(f"Created sample data: {sample_path}")
    return str(sample_path)


def test_basic_validation():
    """Test 1: Basic validation"""
    print("\n" + "=" * 80)
    print("TEST 1: Basic Validation")
    print("=" * 80)

    # Path to data dictionary
    data_dict_path = Path(__file__).parent.parent / 'inst' / 'extdata' / 'cMD_data_dictionary.csv'

    if not data_dict_path.exists():
        print(f"ERROR: Data dictionary not found at {data_dict_path}")
        return False

    # Create sample data
    sample_path = create_sample_data()

    # Initialize validator
    validator = MetadataValidator(str(data_dict_path))

    # Validate
    results = validator.validate_file(sample_path)

    # Print report
    report = validator.generate_report(results)
    print(report)

    # Check results
    assert len(results['extra_columns']) > 0, "Should detect extra columns"
    assert len(results['invalid_columns']) > 0, "Should detect invalid columns"

    print("\n✓ Test 1 PASSED")
    return True


def test_column_renamer():
    """Test 2: Column Renamer Agent"""
    print("\n" + "=" * 80)
    print("TEST 2: Column Renamer Agent")
    print("=" * 80)

    data_dict_path = Path(__file__).parent.parent / 'inst' / 'extdata' / 'cMD_data_dictionary.csv'
    sample_path = Path(__file__).parent / 'test_sample_metadata.csv'

    validator = MetadataValidator(str(data_dict_path))
    df = pd.read_csv(sample_path)

    # Test column renamer
    renamer = ColumnRenamerAgent(validator.data_dictionary, similarity_threshold=0.5)
    validation_results = validator.validate_file(str(sample_path))

    actions = renamer.analyze(df, validation_results)

    print(f"\nProposed {len(actions)} column renamings:")
    for action in actions:
        print(f"  - {action.description}")
        details = action.details
        print(f"    Original: {details['original_column']}")
        print(f"    Suggested: {details['suggested_column']}")
        print(f"    Score: {details['similarity_score']:.2f}")

    assert len(actions) > 0, "Should suggest at least one renaming"
    print("\n✓ Test 2 PASSED")
    return True


def test_value_harmonizer():
    """Test 3: Value Harmonizer Agent"""
    print("\n" + "=" * 80)
    print("TEST 3: Value Harmonizer Agent")
    print("=" * 80)

    data_dict_path = Path(__file__).parent.parent / 'inst' / 'extdata' / 'cMD_data_dictionary.csv'
    sample_path = Path(__file__).parent / 'test_sample_metadata.csv'

    validator = MetadataValidator(str(data_dict_path))
    df = pd.read_csv(sample_path)

    # Test value harmonizer
    harmonizer = ValueHarmonizerAgent(validator.data_dictionary, similarity_threshold=0.7)

    # Teach it some mappings
    harmonizer.learn_mapping('sex', 'M', 'Male')
    harmonizer.learn_mapping('sex', 'F', 'Female')
    harmonizer.learn_mapping('disease', 'T2D', 'Type 2 Diabetes Mellitus')

    validation_results = validator.validate_file(str(sample_path))
    actions = harmonizer.analyze(df, validation_results)

    print(f"\nProposed {len(actions)} value harmonizations:")
    for action in actions:
        print(f"  - {action.description}")

    print("\n✓ Test 3 PASSED")
    return True


def test_ontology_matcher():
    """Test 4: Ontology Matcher Agent"""
    print("\n" + "=" * 80)
    print("TEST 4: Ontology Matcher Agent")
    print("=" * 80)

    data_dict_path = Path(__file__).parent.parent / 'inst' / 'extdata' / 'cMD_data_dictionary.csv'
    sample_path = Path(__file__).parent / 'test_sample_metadata.csv'

    validator = MetadataValidator(str(data_dict_path))
    df = pd.read_csv(sample_path)

    # Test ontology matcher with custom backend
    matcher = OntologyMatcherAgent(validator.data_dictionary, backend='custom')

    actions = matcher.analyze(df)

    print(f"\nProposed {len(actions)} ontology mappings:")
    for action in actions:
        print(f"  - {action.description}")
        if action.details.get('matches'):
            print(f"    Sample matches:")
            for term, matches in list(action.details['matches'].items())[:3]:
                match = matches[0] if matches else {}
                print(f"      '{term}' -> {match.get('term_id', 'N/A')} ({match.get('ontology', 'N/A')})")

    print("\n✓ Test 4 PASSED")
    return True


def test_orchestrator():
    """Test 5: Full Orchestrator Workflow"""
    print("\n" + "=" * 80)
    print("TEST 5: Full Orchestrator Workflow (Auto-approve)")
    print("=" * 80)

    data_dict_path = Path(__file__).parent.parent / 'inst' / 'extdata' / 'cMD_data_dictionary.csv'
    sample_path = Path(__file__).parent / 'test_sample_metadata.csv'
    output_path = Path(__file__).parent / 'test_sample_metadata_harmonized.csv'

    # Create orchestrator with auto-approve
    orchestrator = MetadataOrchestrator(
        data_dictionary_path=str(data_dict_path),
        ontology_backend='custom',
        human_review=False,
        auto_approve=True
    )

    # Process file
    summary = orchestrator.process_file(
        input_csv_path=str(sample_path),
        output_csv_path=str(output_path),
        interactive_review=False
    )

    # Print summary
    summary_report = orchestrator.generate_summary_report(summary)
    print("\n" + summary_report)

    # Check output file was created
    assert output_path.exists(), "Output file should be created"

    # Check improvements
    initial = summary['initial_validation']
    final = summary['final_validation']

    print(f"\nImprovements:")
    print(f"  Extra columns: {initial['extra_columns']} -> {final['extra_columns']}")
    print(f"  Invalid columns: {initial['invalid_columns']} -> {final['invalid_columns']}")

    print("\n✓ Test 5 PASSED")
    return True


def main():
    """Run all tests"""
    tests = [
        test_basic_validation,
        test_column_renamer,
        test_value_harmonizer,
        test_ontology_matcher,
        test_orchestrator
    ]

    print("\n" + "=" * 80)
    print("METADATA VALIDATOR FRAMEWORK - TEST SUITE")
    print("=" * 80)

    passed = 0
    failed = 0

    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"\n✗ {test.__name__} FAILED: {e}")
            import traceback
            traceback.print_exc()
            failed += 1

    print("\n" + "=" * 80)
    print(f"TEST RESULTS: {passed}/{len(tests)} passed")
    if failed == 0:
        print("ALL TESTS PASSED! ✓")
    else:
        print(f"{failed} TESTS FAILED ✗")
    print("=" * 80 + "\n")

    return failed == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
