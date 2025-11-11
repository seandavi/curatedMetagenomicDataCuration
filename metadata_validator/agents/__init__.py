"""
Agentic system for metadata harmonization and correction.

This module provides specialized agents that work together to:
1. Suggest column name corrections
2. Harmonize values to match data dictionary
3. Match terms to ontology concepts
4. Facilitate human review of proposed changes
"""

from .base_agent import BaseAgent, AgentAction
from .column_renamer_agent import ColumnRenamerAgent
from .value_harmonizer_agent import ValueHarmonizerAgent
from .ontology_matcher_agent import OntologyMatcherAgent
from .human_review_agent import HumanReviewAgent

__all__ = [
    'BaseAgent',
    'AgentAction',
    'ColumnRenamerAgent',
    'ValueHarmonizerAgent',
    'OntologyMatcherAgent',
    'HumanReviewAgent'
]
