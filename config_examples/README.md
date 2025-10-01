# PIGSTI Configuration Examples

This directory contains example configuration files to help you get started with PIGSTI.

## Files

- `config.yaml.example` - Main pipeline configuration template
- `samples.tsv.example` - Sample information template
- `config_hops2.0.txt.example` - HOPS configuration template

## Quick Setup

1. Copy the example files to your config directory:
```bash
cp config_examples/*.example config/
```

2. Rename them (remove `.example`):
```bash
cd config/
mv config.yaml.example config.yaml
mv samples.tsv.example samples.tsv
mv config_hops2.0.txt.example config_hops2.0.txt
```

3. Edit the files with your specific paths and parameters

## Configuration Guide

### config.yaml
Main pipeline configuration file containing:
- Database paths (KrakenUniq, BWA indices)
- Analysis parameters
- Quality thresholds
- Resource allocation

### samples.tsv
Tab-separated file with sample information:
- Sample names
- FASTQ file paths (R1 and R2)
- Read group information
- Sequencing run details

### config_hops2.0.txt
HOPS-specific configuration:
- Database paths
- Analysis parameters
- Output settings

## Need Help?

Check the main [README.md](../README.md) for detailed configuration instructions and examples.
