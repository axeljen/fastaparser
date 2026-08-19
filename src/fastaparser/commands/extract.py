"""Extract sequences from FASTA files based on various criteria."""

import sys
import re
import click
import pysam
from Bio.Seq import Seq


def parse_region_string(region_str):
    """Parse a region string in format chr:start-end.
    
    Args:
        region_str: String in format 'chr:start-end' (1-based, inclusive)
        
    Returns:
        Tuple of (chrom, start_0based, end_0based)
    """
    match = re.match(r'^(.+):(\d+)-(\d+)$', region_str)
    if not match:
        raise click.ClickException(
            f"Invalid region format: {region_str}. Expected format: chr:start-end"
        )
    
    chrom = match.group(1)
    start = int(match.group(2)) - 1  # Convert to 0-based
    end = int(match.group(3))  # End is already exclusive in 0-based
    
    if start >= end:
        raise click.ClickException(f"Invalid region: start must be less than end")
    
    return chrom, start, end


def parse_gff(gff_file, feature_id=None):
    """Parse GFF file and extract regions.
    
    Args:
        gff_file: Path to GFF file
        feature_id: Optional feature ID to filter by
        
    Returns:
        Tuple of (list of regions, feature_label)
        Each region is (chrom, start_0based, end_0based, strand)
    """
    regions = []
    feature_label = ""
    
    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom = fields[0]
            feature_type = fields[2]
            start = int(fields[3]) - 1  # GFF is 1-based, convert to 0-based
            end = int(fields[4])  # GFF end is inclusive, our 0-based is exclusive (same value)
            strand = fields[6] if fields[6] in ['+', '-'] else '+'
            attributes = fields[8]
            
            # If feature_id is specified, check if this feature matches
            if feature_id:
                # Check various ID formats in GFF attributes
                if (f"ID={feature_id}" in attributes or 
                    f"Name={feature_id}" in attributes or
                    f"gene_id={feature_id}" in attributes or
                    f"transcript_id={feature_id}" in attributes or
                    feature_id in attributes):
                    regions.append((chrom, start, end, strand))
                    if not feature_label:
                        feature_label = feature_id
            else:
                # No filter, extract all features
                regions.append((chrom, start, end, strand))
    
    if not feature_label:
        feature_label = feature_id if feature_id else "all_features"
    
    return regions, feature_label


def parse_bed(bed_file):
    """Parse BED file and extract regions.
    
    Args:
        bed_file: Path to BED file
        
    Returns:
        List of regions, each is (chrom, start_0based, end_0based, strand)
        BED is already 0-based, half-open
    """
    regions = []
    
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            
            chrom = fields[0]
            start = int(fields[1])  # BED is 0-based
            end = int(fields[2])  # BED end is exclusive (0-based)
            
            # Strand is in 6th column if present
            strand = fields[5] if len(fields) > 5 and fields[5] in ['+', '-'] else '+'
            
            regions.append((chrom, start, end, strand))
    
    return regions


def extract_and_concatenate(fasta, regions, respect_strand):
    """Extract sequences from regions and concatenate them.
    
    Args:
        fasta: pysam.FastaFile object
        regions: List of (chrom, start_0based, end_0based, strand) tuples
        respect_strand: Whether to reverse complement negative strand sequences
        
    Returns:
        Concatenated sequence string
    """
    sequences = []
    
    for chrom, start, end, strand in regions:
        try:
            # Extract sequence (pysam uses 0-based coordinates)
            seq = fasta.fetch(chrom, start, end)
            
            # Reverse complement if on negative strand
            if respect_strand and strand == '-':
                seq = str(Seq(seq).reverse_complement())
            
            sequences.append(seq)
            
        except KeyError:
            raise click.ClickException(
                f"Chromosome/sequence '{chrom}' not found in FASTA file"
            )
        except Exception as e:
            raise click.ClickException(
                f"Error extracting region {chrom}:{start}-{end}: {str(e)}"
            )
    
    return ''.join(sequences)


@click.command()
@click.option(
    "-i",
    "--input",
    "input_file",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True),
    help="Input FASTA file",
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
    "-gff",
    "--gff-file",
    "gff_file",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    default=None,
    help="Optional GFF annotation file",
)
@click.option(
    "-bed",
    "--bed-file",
    "bed_file",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    default=None,
    help="Optional BED file for region extraction",
)
@click.option(
    "-feat",
    "--feature-id",
    "feature_id",
    type=str,
    default=None,
    help="Feature ID to extract (requires --gff-file)",
)
@click.option(
    "-r",
    "--region",
    "region",
    type=str,
    default=None,
    help="Region to extract in format chr:start-end (1-based, inclusive)",
)
@click.option(
    "--strand",
    type=click.Choice(["+", "-"], case_sensitive=True),
    default=None,
    help="Strand orientation (use with --region or --bed-file)",
)
@click.option(
    "--ignore-strand",
    is_flag=True,
    default=False,
    help="Ignore strand information from GFF file (no reverse complement)",
)
@click.option(
    "--header",
    type=str,
    default="subsequence",
    help="Header prefix for output sequence (default: 'subsequence')",
)
def extract(input_file, output_file, gff_file, bed_file, feature_id, region, strand, ignore_strand, header):
    """Extract sequences from a FASTA file.
    
    This command allows you to extract sequences from a FASTA file based on
    various criteria such as GFF annotations, BED regions, or direct coordinates.
    
    Examples:
    
        # Extract a specific region
        fastaparser extract -i input.fa -r chr1:1000-2000
        
        # Extract with reverse complement
        fastaparser extract -i input.fa -r chr1:1000-2000 --strand -
        
        # Extract sequences matching a feature in a GFF file
        fastaparser extract -i input.fa -gff annotations.gff -feat gene123
        
        # Extract sequences from BED regions
        fastaparser extract -i input.fa -bed regions.bed -o output.fa
        
        # Custom header
        fastaparser extract -i input.fa -r chr1:1000-2000 --header my_gene
    """
    # Validate arguments
    if feature_id and not gff_file:
        raise click.UsageError("--feature-id requires --gff-file to be specified")
    
    if strand and not (region or bed_file):
        raise click.UsageError("--strand requires --region or --bed-file to be specified")
    
    if ignore_strand and not gff_file:
        raise click.UsageError("--ignore-strand requires --gff-file to be specified")
    
    # Check that at least one extraction method is specified
    if not any([gff_file, bed_file, region]):
        raise click.UsageError(
            "Please specify at least one of: --gff-file, --bed-file, or --region"
        )
    
    # Determine output stream
    output_stream = sys.stdout if output_file is None else open(output_file, 'w')
    
    try:
        # Open FASTA file with pysam
        fasta = pysam.FastaFile(input_file)
        
        # Collect regions to extract
        regions = []
        feature_label = ""
        respect_strand = False
        
        if gff_file:
            # Parse GFF file and extract regions
            regions, feature_label = parse_gff(gff_file, feature_id)
            if not regions:
                raise click.ClickException(
                    f"No matching features found in GFF file{' for feature ID: ' + feature_id if feature_id else ''}"
                )
            # Use strand from GFF unless ignore_strand is set
            respect_strand = not ignore_strand
            
        elif bed_file:
            # Parse BED file
            regions = parse_bed(bed_file)
            feature_label = f"bed:{bed_file}"
            # If strand is specified, override all regions with that strand
            if strand is not None:
                regions = [(chrom, start, end, strand) for chrom, start, end, _ in regions]
                respect_strand = True
            
        elif region:
            # Parse single region
            chrom, start, end = parse_region_string(region)
            regions = [(chrom, start, end, strand if strand else '+')]
            feature_label = region
            # Only respect strand if explicitly set to '-'
            respect_strand = (strand == '-')
        
        # Extract and concatenate sequences
        concatenated_seq = extract_and_concatenate(fasta, regions, respect_strand)
        
        # Write output
        output_header = f">{header}@{feature_label}\n"
        output_stream.write(output_header)
        
        # Write sequence in 80-character lines
        for i in range(0, len(concatenated_seq), 80):
            output_stream.write(concatenated_seq[i:i+80] + "\n")
        
        if output_file:
            click.echo(f"Extracted {len(concatenated_seq)} bp to {output_file}", err=True)
        else:
            click.echo(f"Extracted {len(concatenated_seq)} bp", err=True)
        
        fasta.close()
        
    except Exception as e:
        raise click.ClickException(str(e))
    finally:
        if output_file and output_stream != sys.stdout:
            output_stream.close()
