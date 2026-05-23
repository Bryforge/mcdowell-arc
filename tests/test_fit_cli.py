import json
from pathlib import Path

import pytest

from mcdowell_arc.cli import main
from mcdowell_arc.fit import load_observations_csv

ROOT = Path(__file__).resolve().parents[1]
SAMPLE = ROOT / "examples" / "sample_webcast.csv"


def test_load_observations_csv_validates_sample():
    obs = load_observations_csv(SAMPLE)
    assert len(obs.t_s) == 11
    assert obs.t_s[0] == 0.0


def test_cli_doctor_outputs_json(capsys):
    assert main(["doctor"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert "mcdowell_arc_version" in payload
    assert "rust_core_available" in payload


def test_cli_fit_smoke(capsys):
    assert main(["fit", str(SAMPLE), "--monte-carlo", "2", "--backend", "python"]) == 0
    payload = json.loads(capsys.readouterr().out)
    assert payload["backend"] == "python"
    assert payload["nominal"]["optimizer_success"] is True
    assert payload["monte_carlo"]["samples_requested"] == 2


def test_bad_csv_missing_required_column(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("t_s,altitude_km\n0,100\n1,101\n2,102\n")
    with pytest.raises(ValueError, match="Missing required CSV columns"):
        load_observations_csv(path)
