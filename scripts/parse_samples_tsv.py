#!/usr/bin/env python3

import csv
import sys
import yaml

def tsv_to_dict(tsv_file, output_yaml):
    samples_dict = {}
    with open(tsv_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sample = row['sample']
            samples_dict[sample] = {
                'r1': row['r1'],
                'r2': row['r2']
            }

    with open(output_yaml, 'w') as out:
        yaml.dump(samples_dict, out, default_flow_style=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python parse_samples_tsv.py samples.tsv samples.yaml")
        sys.exit(1)

    tsv_to_dict(sys.argv[1], sys.argv[2])
