"""
Ontology Matcher Agent

This agent matches terms in the metadata to ontology concepts using
various backends (local OBO files, BioPortal API, or custom matchers).
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Tuple
from difflib import SequenceMatcher
import re
from .base_agent import BaseAgent, AgentAction


class OntologyMatcherAgent(BaseAgent):
    """
    Agent that matches metadata terms to ontology concepts.

    Supports multiple backends:
    - 'pronto': Use pronto library for local OBO ontology files
    - 'bioportal': Use BioPortal REST API
    - 'custom': Use a custom dictionary of term mappings
    """

    def __init__(
        self,
        data_dictionary: Dict[str, Any],
        backend: str = 'custom',
        ontology_config: Optional[Dict[str, Any]] = None
    ):
        """
        Initialize the Ontology Matcher Agent.

        Args:
            data_dictionary: Dictionary mapping column names to DataDictionaryEntry objects
            backend: Backend to use ('pronto', 'bioportal', 'custom')
            ontology_config: Configuration for the backend, e.g.:
                - For 'pronto': {'ontology_path': '/path/to/ontology.obo'}
                - For 'bioportal': {'api_key': 'your_api_key', 'ontologies': ['NCIT', 'DOID']}
                - For 'custom': {'mappings': {...}}
        """
        super().__init__("OntologyMatcherAgent", data_dictionary)
        self.backend = backend
        self.ontology_config = ontology_config or {}
        self.term_cache: Dict[str, List[Dict[str, Any]]] = {}

        # Initialize backend
        if backend == 'pronto':
            self._init_pronto_backend()
        elif backend == 'bioportal':
            self._init_bioportal_backend()
        elif backend == 'custom':
            self._init_custom_backend()
        else:
            raise ValueError(f"Unknown backend: {backend}")

    def _init_pronto_backend(self):
        """Initialize pronto backend for local OBO ontology files"""
        try:
            import pronto
            self.pronto = pronto

            ontology_path = self.ontology_config.get('ontology_path')
            if ontology_path:
                self.ontology = pronto.Ontology(ontology_path)
            else:
                self.ontology = None
                print("Warning: No ontology_path provided for pronto backend")
        except ImportError:
            print("Warning: pronto library not installed. Install with: pip install pronto")
            self.ontology = None

    def _init_bioportal_backend(self):
        """Initialize BioPortal API backend"""
        self.api_key = self.ontology_config.get('api_key')
        self.ontologies = self.ontology_config.get('ontologies', ['NCIT', 'DOID', 'HP'])
        self.base_url = 'https://data.bioportal.org'

        if not self.api_key:
            print("Warning: No API key provided for BioPortal. Set in ontology_config['api_key']")

    def _init_custom_backend(self):
        """Initialize custom mapping backend"""
        # Common disease/condition mappings to standard ontologies
        self.custom_mappings = self.ontology_config.get('mappings', {
            # Disease mappings (simplified examples)
            'Type 2 Diabetes Mellitus': {
                'ontology': 'NCIT',
                'term_id': 'NCIT:C26747',
                'term_label': 'Type 2 Diabetes Mellitus',
                'definition': 'A chronic condition characterized by insulin resistance and relative insulin deficiency.'
            },
            'Type 1 Diabetes Mellitus': {
                'ontology': 'NCIT',
                'term_id': 'NCIT:C2986',
                'term_label': 'Type 1 Diabetes Mellitus',
                'definition': 'An autoimmune disease characterized by destruction of pancreatic beta cells.'
            },
            'Colorectal Carcinoma': {
                'ontology': 'NCIT',
                'term_id': 'NCIT:C2955',
                'term_label': 'Colorectal Carcinoma',
                'definition': 'A malignant epithelial neoplasm arising from the colon or rectum.'
            },
            "Crohn Disease": {
                'ontology': 'DOID',
                'term_id': 'DOID:8778',
                'term_label': "Crohn's disease",
                'definition': 'An inflammatory bowel disease with transmural inflammation.'
            },
            'Rheumatoid Arthritis': {
                'ontology': 'DOID',
                'term_id': 'DOID:7148',
                'term_label': 'Rheumatoid arthritis',
                'definition': 'An autoimmune disease affecting joints.'
            },
            'Healthy': {
                'ontology': 'NCIT',
                'term_id': 'NCIT:C115935',
                'term_label': 'Healthy',
                'definition': 'Being free from illness or disease.'
            },
            # Country mappings
            'United States': {
                'ontology': 'GAZ',
                'term_id': 'GAZ:00002459',
                'term_label': 'United States of America',
                'definition': 'A country in North America.'
            },
            'China': {
                'ontology': 'GAZ',
                'term_id': 'GAZ:00002845',
                'term_label': "People's Republic of China",
                'definition': 'A country in East Asia.'
            },
            # Sex mappings
            'Male': {
                'ontology': 'PATO',
                'term_id': 'PATO:0000384',
                'term_label': 'male',
                'definition': 'An attribute of an individual who has male reproductive organs.'
            },
            'Female': {
                'ontology': 'PATO',
                'term_id': 'PATO:0000383',
                'term_label': 'female',
                'definition': 'An attribute of an individual who has female reproductive organs.'
            },
            # Body site mappings
            'stool': {
                'ontology': 'UBERON',
                'term_id': 'UBERON:0001988',
                'term_label': 'feces',
                'definition': 'Portion of semisolid bodily waste discharged through the anus.'
            },
            'skin': {
                'ontology': 'UBERON',
                'term_id': 'UBERON:0002097',
                'term_label': 'skin of body',
                'definition': 'The organ covering the body that consists of the dermis and epidermis.'
            }
        })

    def analyze(self, df: pd.DataFrame, validation_results: Optional[Dict] = None) -> List[AgentAction]:
        """
        Analyze values and match them to ontology terms.

        Args:
            df: DataFrame to analyze
            validation_results: Optional validation results

        Returns:
            List of proposed ontology matching actions
        """
        proposed_actions = []

        # Target specific columns that should have ontology mappings
        ontology_columns = ['disease', 'country', 'sex', 'body_site', 'target_condition']

        for column in ontology_columns:
            if column not in df.columns:
                continue

            # Get unique values
            unique_values = df[column].dropna().unique()

            matches = {}
            for value in unique_values:
                value_str = str(value)
                if value_str in ['NA', '', 'nan']:
                    continue

                # Handle multiple values (separated by semicolon)
                if ';' in value_str:
                    sub_values = [v.strip() for v in value_str.split(';')]
                    sub_matches = []
                    for sub_val in sub_values:
                        match = self._find_ontology_match(sub_val)
                        if match:
                            sub_matches.append(match)
                    if sub_matches:
                        matches[value_str] = sub_matches
                else:
                    match = self._find_ontology_match(value_str)
                    if match:
                        matches[value_str] = [match]

            if matches:
                action = self.create_action(
                    action_type='add_ontology_mappings',
                    description=f"Add ontology mappings for {len(matches)} terms in '{column}'",
                    details={
                        'column': column,
                        'matches': matches,
                        'new_column_name': f"{column}_ontology"
                    }
                )
                proposed_actions.append(action)

        return proposed_actions

    def _find_ontology_match(self, term: str) -> Optional[Dict[str, Any]]:
        """
        Find ontology match for a term.

        Args:
            term: Term to match

        Returns:
            Dictionary with ontology match information, or None
        """
        # Check cache
        if term in self.term_cache:
            results = self.term_cache[term]
            return results[0] if results else None

        # Search based on backend
        if self.backend == 'pronto':
            result = self._search_pronto(term)
        elif self.backend == 'bioportal':
            result = self._search_bioportal(term)
        elif self.backend == 'custom':
            result = self._search_custom(term)
        else:
            result = None

        # Cache result
        if result:
            self.term_cache[term] = [result]
        else:
            self.term_cache[term] = []

        return result

    def _search_pronto(self, term: str) -> Optional[Dict[str, Any]]:
        """Search for term using pronto backend"""
        if not self.ontology:
            return None

        try:
            # Search for exact match or similar terms
            for ont_term in self.ontology.terms():
                if ont_term.name.lower() == term.lower():
                    return {
                        'ontology': self.ontology.metadata.ontology if hasattr(self.ontology.metadata, 'ontology') else 'OBO',
                        'term_id': ont_term.id,
                        'term_label': ont_term.name,
                        'definition': ont_term.definition if ont_term.definition else '',
                        'synonyms': [str(syn) for syn in ont_term.synonyms] if ont_term.synonyms else []
                    }

            # Fuzzy match
            best_match = None
            best_score = 0
            for ont_term in self.ontology.terms():
                score = SequenceMatcher(None, term.lower(), ont_term.name.lower()).ratio()
                if score > best_score and score > 0.8:
                    best_score = score
                    best_match = ont_term

            if best_match:
                return {
                    'ontology': 'OBO',
                    'term_id': best_match.id,
                    'term_label': best_match.name,
                    'definition': best_match.definition if best_match.definition else '',
                    'similarity_score': best_score
                }

        except Exception as e:
            print(f"Error searching pronto ontology: {e}")

        return None

    def _search_bioportal(self, term: str) -> Optional[Dict[str, Any]]:
        """Search for term using BioPortal API"""
        if not self.api_key:
            return None

        try:
            import requests

            # Use BioPortal search API
            url = f"{self.base_url}/search"
            params = {
                'q': term,
                'ontologies': ','.join(self.ontologies),
                'require_exact_match': 'false',
                'suggest': 'true',
                'pagesize': 1
            }
            headers = {
                'Authorization': f'apikey token={self.api_key}'
            }

            response = requests.get(url, params=params, headers=headers)

            if response.status_code == 200:
                data = response.json()
                if data.get('collection'):
                    result = data['collection'][0]
                    return {
                        'ontology': result.get('links', {}).get('ontology', 'Unknown'),
                        'term_id': result.get('id', ''),
                        'term_label': result.get('prefLabel', ''),
                        'definition': result.get('definition', [''])[0] if result.get('definition') else '',
                        'synonyms': result.get('synonym', [])
                    }

        except Exception as e:
            print(f"Error searching BioPortal: {e}")

        return None

    def _search_custom(self, term: str) -> Optional[Dict[str, Any]]:
        """Search for term using custom mappings"""
        # Exact match
        if term in self.custom_mappings:
            return self.custom_mappings[term]

        # Case-insensitive match
        for key, value in self.custom_mappings.items():
            if key.lower() == term.lower():
                return value

        # Fuzzy match
        best_match = None
        best_score = 0
        for key, value in self.custom_mappings.items():
            score = SequenceMatcher(None, term.lower(), key.lower()).ratio()
            if score > best_score and score > 0.85:
                best_score = score
                best_match = value

        if best_match:
            best_match = best_match.copy()
            best_match['similarity_score'] = best_score
            return best_match

        return None

    def apply_actions(self, df: pd.DataFrame, actions: List[AgentAction]) -> pd.DataFrame:
        """
        Apply approved ontology matching actions.

        Args:
            df: DataFrame to modify
            actions: List of approved actions

        Returns:
            Modified DataFrame with new ontology columns
        """
        df_copy = df.copy()

        for action in actions:
            if action.action_type == 'add_ontology_mappings' and action.approved:
                column = action.details['column']
                matches = action.details['matches']
                new_column = action.details['new_column_name']

                if column in df_copy.columns:
                    # Create ontology column
                    ontology_values = []

                    for value in df_copy[column]:
                        value_str = str(value)
                        if value_str in matches:
                            # Format ontology information
                            match_list = matches[value_str]
                            ontology_ids = [m['term_id'] for m in match_list]
                            ontology_values.append(';'.join(ontology_ids))
                        else:
                            ontology_values.append('')

                    df_copy[new_column] = ontology_values

                    # Also add a label column
                    label_column = f"{column}_ontology_label"
                    label_values = []

                    for value in df_copy[column]:
                        value_str = str(value)
                        if value_str in matches:
                            match_list = matches[value_str]
                            labels = [f"{m['term_label']} ({m['term_id']})" for m in match_list]
                            label_values.append(';'.join(labels))
                        else:
                            label_values.append('')

                    df_copy[label_column] = label_values

                    action.applied = True

        return df_copy

    def add_custom_mapping(self, term: str, ontology: str, term_id: str, term_label: str, definition: str = ''):
        """
        Add a custom term mapping.

        Args:
            term: Original term
            ontology: Ontology abbreviation (e.g., 'NCIT', 'DOID')
            term_id: Ontology term ID
            term_label: Ontology term label
            definition: Term definition
        """
        if self.backend == 'custom':
            self.custom_mappings[term] = {
                'ontology': ontology,
                'term_id': term_id,
                'term_label': term_label,
                'definition': definition
            }

    def export_mappings(self, filepath: str):
        """Export custom mappings to JSON file"""
        import json
        with open(filepath, 'w') as f:
            json.dump(self.custom_mappings, f, indent=2)

    def import_mappings(self, filepath: str):
        """Import custom mappings from JSON file"""
        import json
        with open(filepath, 'r') as f:
            self.custom_mappings.update(json.load(f))
