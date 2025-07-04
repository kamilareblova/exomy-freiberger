import pandas as pd
import argparse
import plotly.graph_objs as go

# Per-base/indel data fields
base_fields = {
    'base': str,
    'count': int,
    'avg_mapping_quality': float,
    'avg_basequality': float,
    'avg_se_mapping_quality': float,
    'num_plus_strand': int,
    'num_minus_strand': int,
    'avg_pos_as_fraction': float,
    'avg_num_mismatches_as_fraction': float,
    'avg_sum_mismatch_qualities': float,
    'num_q2_containing_reads': int,
    'avg_distance_to_q2_start_in_q2_reads': float,
    'avg_clipped_length': float,
    'avg_distance_to_effective_3p_end': float
}

base_map = {
    'A': 5,
    'C': 6,
    'G': 7,
    'T': 8
}

TP_bed = [
    [7572730, 7573106],
    [7573842,	7574141],
    [7576758,	7577244],
    [7577301,	7577675],
    [7578100,	7578627],
    [7579200,	7580022]
]


def count_variants(vcf_file):
    variants = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split()
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4].split(',')[0]
            if len(ref) == 1 and len(alt) == 1 and not chrom.startswith('G') and not chrom.startswith('K') and not chrom.startswith('M'):  # SNPs only
                depth = int(fields[7].split(';')[0].strip('DP='))
                variant = f"{ref}>{alt}"
                variants.append((chrom, pos, variant, depth))
    return variants


def parse_readcounts(readcounts_path):
    # Open the bam-readcount output file and read it line by line
    # Note that the output has no header, so we consume every line
    variants = []
    with open(readcounts_path) as in_fh:
        for line in in_fh:
            # Strip newline from end of line
            line = line.strip()
            # Fields are tab-separated, so split into a list on \t
            fields = line.split('\t')
            # The first four fields contain overall information about the position
            chrom = fields[0]               # Chromosome/reference name
            position = int(fields[1])       # Position (1-based)
            reference_base = fields[2]      # Reference base
            depth = int(fields[3])          # Depth of coverage
            reference_coverage = int(
                fields[base_map[reference_base]].split(':')[1])
            base_ratios = {}
            # The remaining fields are data for each base or indel
            # Iterate over each base/indel
            for base_data_string in fields[4:]:
                # We will store per-base/indel data in a dict
                base_data = {}
                max_base_coverage = reference_coverage

                # Split the base/indel data on ':'
                base_values = base_data_string.split(':')
                # Iterate over each field of the base/indel data
                # Note that this relies on Python 3.5+ to keep the keys in order
                for i, base_field in enumerate(base_fields.keys()):
                    # Store each field of data, converting to the appropriate
                    # data type
                    base_data[base_field] = base_fields[base_field](
                        base_values[i])
                if base_data['base'] == "=" or len(base_data['base']) != 1 or base_data['base'] == "N" or base_data['base'] == reference_base:
                    continue
                if max_base_coverage < base_data['count']:
                    max_base_coverage = base_data['count']
                #base_ratios[f"{base_data['base']}"] = base_data['count']/reference_coverage
                base_ratios[f"{base_data['base']}"] = base_data['count'] / \
                    max_base_coverage

            # append variants
            variants.append({
                'chr': chrom,
                'position': position,
                'ref': reference_base,
                'depth': depth,
                **base_ratios  # unnest object
            })
    return pd.DataFrame(variants)


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(
        description='Plot interactive .html bar plots of bam-readcount output')
    parser.add_argument('--readcount',
                        help='Input file from bam-readcount program')
    parser.add_argument('--name', required=True,
                        help='Name of the output plot file (without extension)')
    args = parser.parse_args()
    # readcounts_path = 'readcounts.txt'
    # readcounts = parse_readcounts(readcounts_path)
    # Count variants for each input file
    readcounts = parse_readcounts(args.readcount)
    for locs in TP_bed:
        df = readcounts[(readcounts['position'] > locs[0]) &
                        (readcounts['position'] < locs[1])]
        traces = []
        for base in ['A', 'C', 'G', 'T']:
            trace = go.Bar(
                x=df['position'],
                y=df[base],
                name=base,
                customdata=df[['position', 'depth']].values
            )
            traces.append(trace)
        # Layout for the plot
        layout = go.Layout(
            title='Base Counts by Position',
            xaxis=dict(
                title='Position',
                tickmode='array',
                tickvals=df['position'],
                ticktext=['{}:{}'.format(row['ref'], row['position'])
                          for _, row in df.iterrows()]
            ),
            yaxis=dict(title='Base Count Ratio', range=[
                       0, 0.04]),  # Adjust range as necessary
            hovermode='x unified', barmode='stack'
        )

        # Create the figure
        fig = go.Figure(data=traces, layout=layout)
        fig.write_html(f"{args.name}-{locs[0]}_{locs[1]}.html")

    readcounts.to_csv("readcounts.csv", index=False)


if __name__ == "__main__":
    main()
