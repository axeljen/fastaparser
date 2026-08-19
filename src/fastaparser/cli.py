"""Main CLI entry point for fastaparser."""

import click
from fastaparser.commands.extract import extract
from fastaparser.commands.filter import filter
from fastaparser.commands.aln_stats import aln_stats


@click.group()
@click.version_option()
def cli():
    """fastaparser - A FASTA parsing and manipulation tool.
    
    Use subcommands to perform various operations on FASTA files.
    """
    pass


# Register subcommands
cli.add_command(extract)
cli.add_command(filter)
cli.add_command(aln_stats)


def main():
    """Entry point for the CLI application."""
    cli()


if __name__ == "__main__":
    main()
