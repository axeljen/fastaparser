# fastaparser

A Python CLI tool for parsing and manipulating FASTA files.

## Features

- **Extract**: Extract sequences from FASTA files using coordinates, GFF annotations, or BED regions
- **Filter**: Filter FASTA alignments based on missing data in sequences and sites

## Installation

Using `uv`:

```bash
# Regular installation
uv pip install .

# Install in development mode (editable)
uv pip install -e .

# Install with dev dependencies for development
uv pip install -e ".[dev]"
```

## Usage

### Extract sequences from FASTA

Extract sequences from a FASTA file:

```bash
# Basic extraction
fastaparser extract -i input.fa -o output.fa

# With GFF annotation (extract by feature)
fastaparser extract -i input.fa -gff annotations.gff -feat gene_id -o output.fa

# With BED regions
fastaparser extract -i input.fa -bed regions.bed -o output.fa

# Output to stdout
fastaparser extract -i input.fa
```

#### Options

- `-i, --input`: Input FASTA file (required)
- `-o, --output`: Output FASTA file (default: stdout)
- `-gff, --gff-file`: Optional GFF annotation file
- `-bed, --bed-file`: Optional BED file
- `-feat, --feature-id`: Feature ID to extract when using GFF
- `-r, --region`: Region to extract in format chr:start-end (1-based, inclusive)
- `--strand`: Strand orientation (+/-)
- `--ignore-strand`: Ignore strand information from GFF file
- `--header`: Header prefix for output sequence

### Filter alignments by missing data

Filter FASTA alignments based on proportion of missing data in sequences and/or sites:

```bash
# Remove sequences with more than 50% missing data
fastaparser filter -i alignment.fa --max-missing-seq 0.5 -o filtered.fa

# Remove sites/columns with more than 20% missing data
fastaparser filter -i alignment.fa --max-missing-site 0.2 -o filtered.fa

# Apply both filters
fastaparser filter -i alignment.fa --max-missing-seq 0.5 --max-missing-site 0.2 -o filtered.fa

# Output to stdout
fastaparser filter -i alignment.fa --max-missing-seq 0.3
```

**Note:** Missing data includes gaps (-), N, X, and ?. Filters are applied sequentially: first sequences are filtered, then sites/columns.

#### Options

- `-i, --input`: Input FASTA alignment file (required)
- `-o, --output`: Output FASTA file (default: stdout)
- `--max-missing-seq`: Maximum proportion of missing data for a sequence to be retained (0.0-1.0)
- `--max-missing-site`: Maximum proportion of missing data for a site/column to be retained (0.0-1.0)

## Development

Install in development mode with dev dependencies:

```bash
uv pip install -e ".[dev]"
```

Run tests:

```bash
pytest
```

Format code:

```bash
black src/
ruff check src/
```

## Adding New Subcommands

To add a new subcommand:

1. Create a new file in `src/fastaparser/commands/`
2. Define your command using Click decorators
3. Register it in `src/fastaparser/cli.py`

Example:

```python
# src/fastaparser/commands/mynewcommand.py
import click

@click.command()
@click.option('-i', '--input', required=True, help='Input file')
def mynewcommand(input):
    """Description of my new command."""
    click.echo(f"Processing {input}")

# src/fastaparser/cli.py
from fastaparser.commands.mynewcommand import mynewcommand
cli.add_command(mynewcommand)
```
