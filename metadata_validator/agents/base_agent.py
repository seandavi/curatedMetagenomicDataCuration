"""
Base Agent class for the metadata harmonization system.

All specialized agents inherit from this base class.
"""

import pandas as pd
from abc import ABC, abstractmethod
from typing import Dict, List, Any, Optional
from dataclasses import dataclass, field
from datetime import datetime
import json


@dataclass
class AgentAction:
    """Represents a single action taken by an agent"""
    agent_name: str
    action_type: str
    timestamp: str
    description: str
    details: Dict[str, Any]
    approved: bool = False
    applied: bool = False

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization"""
        return {
            'agent_name': self.agent_name,
            'action_type': self.action_type,
            'timestamp': self.timestamp,
            'description': self.description,
            'details': self.details,
            'approved': self.approved,
            'applied': self.applied
        }


class BaseAgent(ABC):
    """
    Base class for all metadata harmonization agents.

    Provides common functionality for:
    - Tracking actions and history
    - Logging changes
    - Interacting with data dictionary
    """

    def __init__(self, name: str, data_dictionary: Dict[str, Any]):
        """
        Initialize the agent.

        Args:
            name: Name of the agent
            data_dictionary: Dictionary mapping column names to DataDictionaryEntry objects
        """
        self.name = name
        self.data_dictionary = data_dictionary
        self.actions: List[AgentAction] = []
        self.enabled = True

    @abstractmethod
    def analyze(self, df: pd.DataFrame, validation_results: Optional[Dict] = None) -> List[AgentAction]:
        """
        Analyze the dataframe and propose actions.

        Args:
            df: DataFrame to analyze
            validation_results: Optional validation results from MetadataValidator

        Returns:
            List of proposed AgentAction objects
        """
        pass

    @abstractmethod
    def apply_actions(self, df: pd.DataFrame, actions: List[AgentAction]) -> pd.DataFrame:
        """
        Apply approved actions to the dataframe.

        Args:
            df: DataFrame to modify
            actions: List of approved actions to apply

        Returns:
            Modified DataFrame
        """
        pass

    def create_action(self, action_type: str, description: str, details: Dict[str, Any]) -> AgentAction:
        """
        Create a new AgentAction.

        Args:
            action_type: Type of action (e.g., 'rename_column', 'harmonize_value')
            description: Human-readable description
            details: Dictionary with action-specific details

        Returns:
            New AgentAction object
        """
        action = AgentAction(
            agent_name=self.name,
            action_type=action_type,
            timestamp=datetime.now().isoformat(),
            description=description,
            details=details
        )
        self.actions.append(action)
        return action

    def get_action_history(self) -> List[Dict[str, Any]]:
        """Get history of all actions taken by this agent"""
        return [action.to_dict() for action in self.actions]

    def get_approved_actions(self) -> List[AgentAction]:
        """Get all approved but not yet applied actions"""
        return [action for action in self.actions if action.approved and not action.applied]

    def clear_history(self):
        """Clear the action history"""
        self.actions.clear()

    def save_history(self, filepath: str):
        """Save action history to JSON file"""
        with open(filepath, 'w') as f:
            json.dump([action.to_dict() for action in self.actions], f, indent=2)

    def load_history(self, filepath: str):
        """Load action history from JSON file"""
        with open(filepath, 'r') as f:
            data = json.load(f)
            self.actions = [
                AgentAction(
                    agent_name=a['agent_name'],
                    action_type=a['action_type'],
                    timestamp=a['timestamp'],
                    description=a['description'],
                    details=a['details'],
                    approved=a.get('approved', False),
                    applied=a.get('applied', False)
                )
                for a in data
            ]
