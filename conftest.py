"""pytest configuration: add project root to module search path."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
