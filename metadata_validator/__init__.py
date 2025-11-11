"""
Metadata Validator Framework for Curated Metagenomic Data

A comprehensive Python framework for validating, harmonizing, and enriching
sample metadata against a data dictionary with agentic assistance.

Main components:
- MetadataValidator: Core validation engine
- MetadataOrchestrator: Workflow coordinator
- Agents: Specialized harmonization agents (column renaming, value harmonization, ontology matching)
"""

from .validator import MetadataValidator, ValidationResult, DataDictionaryEntry
from .orchestrator import MetadataOrchestrator
from .agents import (
    BaseAgent,
    AgentAction,
    ColumnRenamerAgent,
    ValueHarmonizerAgent,
    OntologyMatcherAgent,
    HumanReviewAgent
)

__version__ = '0.1.0'

__all__ = [
    # Core validation
    'MetadataValidator',
    'ValidationResult',
    'DataDictionaryEntry',

    # Orchestration
    'MetadataOrchestrator',

    # Agents
    'BaseAgent',
    'AgentAction',
    'ColumnRenamerAgent',
    'ValueHarmonizerAgent',
    'OntologyMatcherAgent',
    'HumanReviewAgent'
]
