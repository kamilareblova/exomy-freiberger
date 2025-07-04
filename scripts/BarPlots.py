import pandas as pd
import collections
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

order = [
    'A>C', 'A>G', 'A>T',
    'C>A', 'C>G', 'C>T',
    'G>A', 'G>C', 'G>T',
    'T>A', 'T>C', 'T>G']

# Function to count variant types


def count_variants(vcf_file):
    variant_counts = collections.Counter()
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.split()
            chrom = fields[0]
            ref = fields[3]
            alt = fields[4]
            # skip funny chromosomes and indels
            if len(ref) == 1 and len(alt) == 1 and not chrom.startswith('G') and not chrom.startswith('K') and not chrom.startswith('M'):
                try:
                    # parse either of the variantcallers
                    depth = int(fields[11].split(':')[2])
                except:
                    # vardict needs extra care
                    depth = int(fields[9].split(':')[1].split(',')[1])

                variant = "{}>{}".format(ref, alt)
                variant_counts[variant] += depth
    return variant_counts


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Count variant types in VCF files and plot histogram.')
    parser.add_argument('vcf_files', nargs='+', metavar='VCF_FILE',
                        help='List of input VCF files')
    parser.add_argument('--name', required=True,
                        help='Name of the output plot file (without extension)')
    args = parser.parse_args()

    all_variant_counts = collections.Counter()
    for vcf_file in args.vcf_files:
        variant_counts = count_variants(vcf_file)
        all_variant_counts.update(variant_counts)

    variant_df = pd.DataFrame(
        list(all_variant_counts.items()), columns=['Variant', 'Count'])

    variant_df.to_csv(args.name + "_counts.csv", index=False)

    plt.figure(figsize=(12, 6))
    sns.barplot(x='Variant', y='Count', data=variant_df, palette='muted', order=order)
    plt.title('Counts of Variant Types in ' + args.name)
    plt.xlabel('Variant Type')
    plt.ylabel('Count')
    plt.xticks(rotation=45)
    plt.savefig(args.name + '.pdf')


if __name__ == "__main__":
    main()
