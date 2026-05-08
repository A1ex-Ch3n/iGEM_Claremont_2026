"""
pytest conftest: add processes/ to sys.path so rbp_pipeline can be imported.
将 processes/ 目录添加到 sys.path，以便导入 rbp_pipeline。
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
