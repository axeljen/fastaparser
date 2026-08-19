"""Filter sequences from FASTA alignments based on missing data."""

import sys
import click
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def count_missing(sequence):
    """Count missing/ambiguous characters in a sequence.
    
    Missing data includes: gaps (-), N, X, and ?
    
    Args:
        sequence: String or Seq object
        
    Returns:
        Number of missing characters
    """
    seq_str = str(sequence).upper()
    missing_chars = {'-', 'N', 'X', '?'}
    return sum(1 for c in seq_str if c in missing_chars)


def filter_sequences_by_missing(records, max_missing_prop):
    """Filter sequences based on proportion of missing data.
    
    Args:
        records: List of SeqRecord objects
        max_missing_prop: Maximum proportion of missing data (0.0 to 1.0)
        
    Returns:
        Filtered list of SeqRecord objects
    """
    filtered = []
    
    for record in records:
        seq_len = len(record.seq)
        if seq_len == 0:
            continue
            
        missing_count = count_missing(record.seq)
        missing_prop = missing_count / seq_len
        
        if missing_prop <= max_missing_prop:
            filtered.append(record)
    
    return filtered


def filter_sites_by_missing(records, max_missing_prop):
    """Filter alignment columns/sites based on proportion of missing data.
    
    Args:
        records: List of SeqRecord objects (aligned sequences)
        max_missing_prop: Maximum proportion of missing data (0.0 to 1.0)
        
    Returns:
        List of SeqRecord objects with filtered columns
    """
    if not records:
        return records
    
    # Get alignment length
    aln_length = len(records[0].seq)
    
    # Check all sequences have same length (alignment requirement)
    if not all(len(rec.seq) == aln_length for rec in records):
        raise click.ClickException(
            "Input sequences must be aligned (all same length) to filter by site"
        )
    
    num_sequences = len(records)
    
    # Determine which columns to keep
    columns_to_keep = []
    
    for col_idx in range(aln_length):
        # Count missing data in this column
        missing_count = sum(
            1 for rec in records 
            if str(rec.seq[col_idx]).upper() in {'-', 'N', 'X', '?'}
        )
        
        missing_prop = missing_count / num_sequences
        
        if missing_prop <= max_missing_prop:
            columns_to_keep.append(col_idx)
    
    # Build filtered sequences
    filtered_records = []
    for record in records:
        seq_str = str(record.seq)
        filtered_seq = ''.join(seq_str[i] for i in columns_to_keep)
        
        new_record = SeqRecord(
            Seq(filtered_seq),
            id=record.id,
            description=record.description,
            name=record.name
        )
        filtered_records.append(new_record)
    
    return filtered_records


@click.command()
@click.option(
    "-i",
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Input FASTA alignment file",
)
@click.option(
    "-o",
    "--output",
    "output_file",
    type=click.Path(dir_okay=False, writable=True),
    default=None,
    help="Output FASTA file (default: stdout)",
)
@click.option(
    "--max-missing-seq",
    "max_missing_seq",
    type=click.FloatRange(0.0, 1.0),
    default=None,
    help="Maximum proportion of missing data for a sequence to be retained (0.0-1.0)",
)
@click.option(
    "--max-missing-site",
    "max_missing_site",
    type=click.FloatRange(0.0, 1.0),
    default=None,
    help="Maximum proportion of missing data for a site/column to be retained (0.0-1.0)",
)
def filter(input_file, output_file, max_missing_seq, max_missing_site):
    """Filter FASTA alignment based on missing data.
    
    This command filters sequences and/or alignment columns based on the
    proportion of missing data. Missing data includes gaps (-), N, X, and ?.
    
    Filters are applied in order:
    1. First, sequences exceeding --max-missing-seq are removed
    2. Then, sites/columns exceeding --max-missing-site are removed
    
    Examples:
    
        # Remove sequences with more than 50% missing data
        fastaparser filter -i alignment.fa --max-missing-seq 0.5 -o filtered.fa
        
        # Remove sites with more than 20% missing data
        fastaparser filter -i alignment.fa --max-missing-site 0.2 -o filtered.fa
        
        # Apply both filters
        fastaparser filter -i alignment.fa --max-missing-seq 0.5 --max-missing-site 0.2 -o filtered.fa
    """
    # Validate that at least one filter is specified
    if max_missing_seq is None and max_missing_site is None:
        raise click.UsageError(
            "Please specify at least one filter: --max-missing-seq or --max-missing-site"
        )
    
    # Determine output stream
    output_stream = sys.stdout if output_file is None else open(output_file, 'w')
    
    try:
        # Read input sequences
        records = list(SeqIO.parse(input_file, "fasta"))
        
        if not records:
            raise click.ClickException("No sequences found in input file")
        
        initial_seq_count = len(records)
        initial_length = len(records[0].seq) if records else 0
        
        # Apply sequence filter if specified
        if max_missing_seq is not None:
            records = filter_sequences_by_missing(records, max_missing_seq)
            click.echo(
                f"Sequence filter: retained {len(records)}/{initial_seq_count} sequences",
                err=True
            )
        
        if not records:
            raise click.ClickException(
                "No sequences remain after filtering. Try relaxing --max-missing-seq threshold."
            )
        
        # Apply site filter if specified
        if max_missing_site is not None:
            records = filter_sites_by_missing(records, max_missing_site)
            final_length = len(records[0].seq) if records else 0
            click.echo(
                f"Site filter: retained {final_length}/{initial_length} sites",
                err=True
            )
        
        if not records or (records and len(records[0].seq) == 0):
            raise click.ClickException(
                "No sites remain after filtering. Try relaxing --max-missing-site threshold."
            )
        
        # Write output
        SeqIO.write(records, output_stream, "fasta")
        
        if output_file:
            click.echo(
                f"Filtered alignment written to {output_file}: "
                f"{len(records)} sequences, {len(records[0].seq)} sites",
                err=True
            )
        else:
            click.echo(
                f"Filtered alignment: {len(records)} sequences, {len(records[0].seq)} sites",
                err=True
            )
        
    except Exception as e:
        if not isinstance(e, click.ClickException):
            raise click.ClickException(str(e))
        raise
    finally:
        if output_file and output_stream != sys.stdout:
            output_stream.close()
