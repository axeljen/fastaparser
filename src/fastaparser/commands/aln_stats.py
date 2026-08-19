"""Get statistics from a FASTA alignment."""

import sys
import click
from Bio import SeqIO


@click.command()
@click.option(
    "-i",
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Input FASTA alignment file",
)
def aln_stats(input_file):
    """Get statistics from a FASTA alignment.
    
    Outputs the number of sequences and total number of bases in the alignment.
    
    Exits with an error if sequences have different lengths (not an alignment).
    
    Examples:
    
        # Get alignment statistics
        fastaparser aln_stats -i alignment.fa
    """
    try:
        # Read input sequences
        records = list(SeqIO.parse(input_file, "fasta"))
        
        if not records:
            raise click.ClickException("No sequences found in input file")
        
        # Get number of sequences
        num_sequences = len(records)
        
        # Check that all sequences have the same length
        seq_lengths = [len(rec.seq) for rec in records]
        
        if len(set(seq_lengths)) > 1:
            # Sequences have different lengths - not an alignment
            length_info = ", ".join(f"{rec.id}:{len(rec.seq)}" for rec in records[:5])
            if len(records) > 5:
                length_info += ", ..."
            raise click.ClickException(
                f"Sequences have different lengths (not an alignment).\n"
                f"First few sequences: {length_info}"
            )
        
        # All sequences have the same length
        aln_length = seq_lengths[0]
        total_bases = num_sequences * aln_length
        
        # Output statistics
        click.echo(f"Number of sequences: {num_sequences}\nAlignment length: {aln_length}\nTotal bases: {total_bases}")

    except Exception as e:
        if not isinstance(e, click.ClickException):
            raise click.ClickException(str(e))
        raise
