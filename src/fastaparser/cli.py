"""Main CLI entry point for fastaparser."""

import click
from fastaparser.commands.extract import extract
from fastaparser.commands.filter import filter


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


def main():
    """Entry point for the CLI application."""
    cli()


if __name__ == "__main__":
    main()
