"""CLI tests: --help snapshot, version, offline benchmark, and an e2e `run`."""

import json

import pytest
from typer.testing import CliRunner

from proteoform_regions import __version__
from proteoform_regions.cli import app
from tests.synthetic import write_cohort

runner = CliRunner()

COMMANDS = ["harmonize", "map", "compute", "background", "stitch", "bias", "run", "benchmark"]


class TestSurface:
    def test_version(self):
        result = runner.invoke(app, ["--version"])
        assert result.exit_code == 0
        assert __version__ in result.output

    def test_help_lists_every_subcommand(self):
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        for command in COMMANDS:
            assert command in result.output

    @pytest.mark.parametrize("command", COMMANDS)
    def test_subcommand_help_snapshot(self, command):
        result = runner.invoke(app, [command, "--help"])
        assert result.exit_code == 0
        assert "--help" in result.output
        # paper defaults are visible in every stage's help
        if command in ("stitch", "run"):
            assert "--gap" in result.output
        if command in ("bias", "run"):
            assert "--seed" in result.output

    def test_unknown_command_fails(self):
        assert runner.invoke(app, ["warp"]).exit_code != 0


class TestBenchmark:
    def test_offline_benchmark(self):
        result = runner.invoke(app, ["benchmark"])
        assert result.exit_code == 0
        assert "benchmark rows: 16" in result.output
        # paper-known figures: poly-E outlier 0.420; median non-outlier agreement
        assert "max|d|=0.4200" in result.output


class TestEndToEnd:
    def test_run_command(self, tmp_path):
        cohort_yaml, data_dir, cache_dir = write_cohort(tmp_path)
        out_dir = tmp_path / "run"
        result = runner.invoke(
            app,
            [
                "run",
                "--study",
                str(cohort_yaml),
                "--out-dir",
                str(out_dir),
                "--data-dir",
                str(data_dir),
                "--uniprot-cache",
                str(cache_dir),
                "--no-download",
                "--offline",
                "--b-perm",
                "50",
                "--b-boot",
                "50",
            ],
        )
        assert result.exit_code == 0, result.output
        assert (out_dir / "manifest.sha256").exists()
        payload = json.loads((out_dir / "processed_inference_results.json").read_text())
        assert payload["k"] == 2

    def test_stage_commands_compose(self, tmp_path):
        """Individual subcommands chain on intermediate CSVs like the notebooks."""
        cohort_yaml, data_dir, cache_dir = write_cohort(tmp_path)
        work = tmp_path / "stages"
        work.mkdir()
        common = ["--evidence", str(work / "evidence.csv")]

        r1 = runner.invoke(
            app,
            [
                "harmonize",
                "--study",
                str(cohort_yaml),
                "--data-dir",
                str(data_dir),
                "--no-download",
                "--out",
                str(work / "evidence.csv"),
            ],
        )
        assert r1.exit_code == 0, r1.output

        r2 = runner.invoke(
            app,
            [
                "map",
                *common,
                "--uniprot-cache",
                str(cache_dir),
                "--offline",
                "--out-map",
                str(work / "map.csv"),
                "--out-properties",
                str(work / "props.csv"),
            ],
        )
        assert r2.exit_code == 0, r2.output

        r3 = runner.invoke(
            app,
            [
                "compute",
                *common,
                "--map",
                str(work / "map.csv"),
                "--proteins",
                str(work / "props.csv"),
                "--out",
                str(work / "features.csv"),
            ],
        )
        assert r3.exit_code == 0, r3.output

        r4 = runner.invoke(
            app,
            [
                "background",
                "--proteins",
                str(work / "props.csv"),
                "--features",
                str(work / "features.csv"),
                "--out",
                str(work / "bg.csv"),
            ],
        )
        # properties CSV written by `map` drops sequences -> explicit, honest failure
        assert r4.exit_code == 2

        r5 = runner.invoke(
            app,
            [
                "background",
                "--uniprot-cache",
                str(cache_dir),
                "--features",
                str(work / "features.csv"),
                "--out",
                str(work / "bg.csv"),
            ],
        )
        assert r5.exit_code == 0, r5.output
        assert (work / "bg.csv").exists()

        r6 = runner.invoke(
            app,
            [
                "stitch",
                "--features",
                str(work / "features.csv"),
                "--map",
                str(work / "map.csv"),
                "--uniprot-cache",
                str(cache_dir),
                "--out",
                str(work / "regions.csv"),
            ],
        )
        assert r6.exit_code == 0, r6.output
        regions = (work / "regions.csv").read_text()
        assert "region_supported" in regions

        r7 = runner.invoke(
            app,
            [
                "bias",
                "--features",
                str(work / "features.csv"),
                "--background",
                str(work / "bg.csv"),
                "--no-rq1",
                "--out",
                str(work / "bias.json"),
            ],
        )
        assert r7.exit_code == 0, r7.output
        payload = json.loads((work / "bias.json").read_text())
        assert payload["k"] == 2
