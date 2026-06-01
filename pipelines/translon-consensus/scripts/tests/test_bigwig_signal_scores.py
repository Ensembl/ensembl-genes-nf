import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPT_DIR))

from transcode_bigwig_signal_scores import ScoreConfig, periodicity_score, score_signal_array


def test_periodicity_score_uses_calibrated_annotated_frame():
    body = np.zeros(30, dtype=float)
    body[np.arange(len(body)) % 3 == 2] = 1.0

    assert periodicity_score(body, start_offset=0) == 0.0
    assert periodicity_score(body, start_offset=1) == 1.0


def test_score_signal_array_applies_psite_offset():
    arr = np.zeros(60, dtype=float)
    body_positions = np.arange(30) % 3 == 2
    arr[15:45][body_positions] = 1.0

    assert score_signal_array(arr, ScoreConfig(psite_offset=0))["periodicity"] == 0.0
    assert score_signal_array(arr, ScoreConfig(psite_offset=1))["periodicity"] == 1.0
