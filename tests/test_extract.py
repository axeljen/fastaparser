"""Tests for the extract command."""

import pytest
from click.testing import CliRunner
from fastaparser.cli import cli


def test_extract_requires_input():
    """Test that extract command requires -i flag."""
    runner = CliRunner()
    result = runner.invoke(cli, ["extract"])
    assert result.exit_code != 0
    assert "Missing option" in result.output or "required" in result.output.lower()


def test_extract_help():
    """Test that extract command shows help."""
    runner = CliRunner()
    result = runner.invoke(cli, ["extract", "--help"])
    assert result.exit_code == 0
    assert "Extract sequences from a FASTA file" in result.output


def test_extract_feature_without_gff():
    """Test that --feature-id requires --gff-file."""
    runner = CliRunner()
    with runner.isolated_filesystem():
        # Create a dummy input file
        with open("test.fa", "w") as f:
            f.write(">seq1\nATCG\n")
        
        result = runner.invoke(cli, ["extract", "-i", "test.fa", "-feat", "gene123"])
        assert result.exit_code != 0
        assert "feature-id requires --gff-file" in result.output
