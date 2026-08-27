#!/usr/bin/env python3
"""
Filter ENA sheets for RANDOM library selection, download reads, and create PIGSTI samples.tsv.

This script:
1. Reads one or more ENA file reports (TSV format)
2. Filters for library_selection == "RANDOM"
3. Downloads FASTQ files (handles both PAIRED and SINGLE end)
4. Creates a PIGSTI-compatible samples.tsv file

Usage:
    python ena_to_pigsti_samplesheet.py <ena_file1> [<ena_file2> ...] --output-dir <output_dir> [--samples-tsv <output_tsv>] [--dry-run]

Example:
    python ena_to_pigsti_samplesheet.py filereport_read_run_PRJEB26011.tsv filereport_read_run_PRJEB31621.tsv --output-dir ./raw_data --samples-tsv config/samples.tsv
"""

import os
import sys
import argparse
import pandas as pd
import subprocess
from pathlib import Path
from urllib.parse import urlparse
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import download functions from existing download_ena_fastq.py
# We'll reuse the download logic but add filtering for RANDOM library_selection

def check_aria2c():
    """Check if aria2c is available."""
    try:
        subprocess.run(['aria2c', '--version'], capture_output=True, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def download_file_aria2c(url, output_path, max_connections=4, max_retries=5, retry_wait=3):
    """Download a file using aria2c (parallel downloader)."""
    cmd = [
        'aria2c',
        '-x', str(max_connections),
        '-s', str(max_connections),
        '--min-split-size', '10M',
        '--max-tries', str(max_retries),
        '--retry-wait', str(retry_wait),
        '--timeout', '60',
        '--connect-timeout', '20',
        '--continue',
        '--file-allocation', 'none',
        '--allow-overwrite',
        '--auto-file-renaming=false',
        '--max-connection-per-server', str(max_connections),
        '--split', str(max_connections),
        '--console-log-level=warn',
        '--summary-interval=0',
        '-d', str(output_path.parent),
        '-o', output_path.name,
        url
    ]
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=10800)
        if output_path.exists() and output_path.stat().st_size > 0:
            return True, None
        else:
            return False, "File downloaded but is empty or missing"
    except subprocess.TimeoutExpired:
        return False, "Download timeout (exceeded 3 hours)"
    except subprocess.CalledProcessError as e:
        error_parts = []
        if e.stderr and e.stderr.strip():
            error_parts.append(f"stderr: {e.stderr.strip()}")
        if e.stdout and e.stdout.strip():
            error_parts.append(f"stdout: {e.stdout.strip()}")
        if not error_parts:
            if not output_path.exists():
                error_parts.append(f"File not created (exit code {e.returncode})")
            else:
                error_parts.append(f"Exit code: {e.returncode}")
        
        error_msg = " | ".join(error_parts) if error_parts else f"Exit code: {e.returncode}"
        return False, f"aria2c failed (code {e.returncode}): {error_msg[:200]}"

def download_file_single(url, output_path, max_retries=5, ftp_url=None):
    """Download a file from URL with retries (single-threaded)."""
    last_error = None
    urls_to_try = [url]
    if ftp_url:
        urls_to_try.append(ftp_url)
    
    for attempt in range(max_retries):
        if attempt > 0:
            delay = min(10 * attempt, 30)
            time.sleep(delay)
        
        for current_url in urls_to_try:
            try:
                if current_url.startswith('ftp://'):
                    http_url = current_url.replace('ftp://', 'http://', 1)
                    cmd = ['wget', '-O', str(output_path), '--continue', '--timeout=120', '--tries=1', '--waitretry=5', http_url]
                elif current_url.startswith('http://') or current_url.startswith('https://'):
                    cmd = ['wget', '-O', str(output_path), '--continue', '--timeout=120', '--tries=1', '--waitretry=5', current_url]
                else:
                    http_url = 'http://' + current_url
                    cmd = ['wget', '-O', str(output_path), '--continue', '--timeout=120', '--tries=1', '--waitretry=5', http_url]
                
                result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=10800)
                if output_path.exists() and output_path.stat().st_size > 0:
                    return True, None
                else:
                    if output_path.exists():
                        output_path.unlink()
                    last_error = "File downloaded but is empty"
                    continue
            except (subprocess.TimeoutExpired, subprocess.CalledProcessError, FileNotFoundError) as e:
                if output_path.exists():
                    try:
                        output_path.unlink()
                    except:
                        pass
                if isinstance(e, FileNotFoundError):
                    # Try curl
                    try:
                        if current_url.startswith('ftp://'):
                            cmd = ['curl', '-o', str(output_path), '-C', '-', '--max-time', '10800', '--retry', '3', '--retry-delay', '5', '--ftp-pasv', current_url]
                        else:
                            http_url = current_url.replace('ftp://', 'http://', 1) if current_url.startswith('ftp://') else current_url
                            cmd = ['curl', '-L', '-o', str(output_path), '-C', '-', '--max-time', '10800', '--retry', '3', '--retry-delay', '5', '--fail', http_url]
                        result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=10800)
                        if output_path.exists() and output_path.stat().st_size > 0:
                            return True, None
                    except:
                        pass
                last_error = f"Download failed: {str(e)[:100]}"
                continue
    
    return False, last_error or "Download failed after all retries"

def download_file_wrapper(url, output_path, use_aria2c=False, ftp_url=None, connections_per_file=2):
    """Wrapper function for downloading a single file."""
    if output_path.exists() and output_path.stat().st_size > 0:
        return (True, f"Already exists: {output_path.name}", None)
    
    if use_aria2c:
        success, error_msg = download_file_aria2c(url, output_path, max_connections=connections_per_file, max_retries=5, retry_wait=5)
        if not success:
            time.sleep(2)
            print(f"      Trying fallback (wget/curl) for {output_path.name}...")
            fallback_success, fallback_error = download_file_single(url, output_path, max_retries=5, ftp_url=ftp_url)
            if fallback_success:
                return (True, f"Downloaded (fallback): {output_path.name}", None)
            else:
                return (False, f"Failed: {output_path.name}", f"aria2c: {error_msg}; fallback: {fallback_error}")
    else:
        success, error_msg = download_file_single(url, output_path, max_retries=5, ftp_url=ftp_url)
        if not success:
            error_msg = error_msg or "Download failed"
    
    if success:
        return (True, f"Downloaded: {output_path.name}", None)
    else:
        return (False, f"Failed: {output_path.name}", error_msg)

def format_bytes(bytes_value):
    """Convert bytes to human-readable format."""
    if pd.isna(bytes_value) or bytes_value == 0:
        return "Unknown"
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if bytes_value < 1024.0:
            return f"{bytes_value:.2f} {unit}"
        bytes_value /= 1024.0
    return f"{bytes_value:.2f} PB"

def parse_ena_file(ena_file, filter_random=True):
    """
    Parse ENA file report and extract download information.
    Filters for library_selection == "RANDOM" if filter_random is True.
    Returns samples_data and total_size_bytes.
    """
    df = pd.read_csv(ena_file, sep='\t')
    
    # Filter for RANDOM library selection
    if filter_random:
        if 'library_selection' not in df.columns:
            print(f"[WARNING] No 'library_selection' column in {ena_file}, skipping filter")
        else:
            initial_count = len(df)
            df = df[df['library_selection'].str.upper() == 'RANDOM'].copy()
            filtered_count = len(df)
            print(f"  Filtered to {filtered_count} samples with library_selection=RANDOM (from {initial_count} total)")
    
    samples_data = []
    total_size_bytes = 0
    
    for idx, row in df.iterrows():
        run_accession = row['run_accession']
        sample_accession = row.get('sample_accession', '')
        study_accession = row.get('study_accession', '')
        library_name = row.get('library_name', '')
        library_layout = row.get('library_layout', '').upper()  # PAIRED or SINGLE
        
        # Get file size from sra_bytes column if available
        sra_bytes = row.get('sra_bytes', 0)
        try:
            if pd.notna(sra_bytes) and sra_bytes != '':
                file_size = int(float(sra_bytes))
            else:
                file_size = 0
        except (ValueError, TypeError):
            file_size = 0
        
        # Check for fastq_ftp first (newer format), then submitted_ftp (older format)
        fastq_ftp = row.get('fastq_ftp', '')
        submitted_ftp = row.get('submitted_ftp', '')
        
        # Prefer fastq_ftp, fallback to submitted_ftp
        ftp_urls_str = fastq_ftp if pd.notna(fastq_ftp) and fastq_ftp != '' else submitted_ftp
        
        # Skip if no FTP URLs
        if pd.isna(ftp_urls_str) or ftp_urls_str == '':
            print(f"[WARNING] No FTP URLs for {run_accession}, skipping...")
            continue
        
        # Parse FTP URLs (semicolon-separated for paired, single URL for single-end)
        ftp_urls = [url.strip() for url in str(ftp_urls_str).split(';') if url.strip()]
        
        # Find R1 and R2 files
        r1_url = None
        r1_url_ftp = None
        r2_url = None
        r2_url_ftp = None
        
        if library_layout == 'SINGLE':
            # Single-end: only one file
            if len(ftp_urls) > 0:
                url = ftp_urls[0].strip()
                original_url = url
                if original_url.startswith('ftp://'):
                    http_url = original_url.replace('ftp://', 'http://', 1)
                    ftp_url = original_url
                elif not original_url.startswith('http'):
                    http_url = 'http://' + original_url
                    ftp_url = 'ftp://' + original_url
                else:
                    http_url = original_url
                    ftp_url = original_url.replace('http://', 'ftp://', 1) if original_url.startswith('http://') else None
                
                r1_url = http_url
                r1_url_ftp = ftp_url
        else:
            # Paired-end: find R1 and R2
            for url in ftp_urls:
                original_url = url.strip()
                if original_url.startswith('ftp://'):
                    http_url = original_url.replace('ftp://', 'http://', 1)
                    ftp_url = original_url
                elif not original_url.startswith('http'):
                    http_url = 'http://' + original_url
                    ftp_url = 'ftp://' + original_url
                else:
                    http_url = original_url
                    ftp_url = original_url.replace('http://', 'ftp://', 1) if original_url.startswith('http://') else None
                
                filename = os.path.basename(original_url)
                if '_R1_' in filename or '_r1_' in filename or filename.endswith('_1.fastq.gz') or filename.endswith('_1.fq.gz'):
                    r1_url = http_url
                    r1_url_ftp = ftp_url
                elif '_R2_' in filename or '_r2_' in filename or filename.endswith('_2.fastq.gz') or filename.endswith('_2.fq.gz'):
                    r2_url = http_url
                    r2_url_ftp = ftp_url
        
        # Add file size to total (sra_bytes represents total size for the run, whether single or paired)
        if file_size > 0:
            total_size_bytes += file_size
        
        if not r1_url:
            print(f"[WARNING] No R1 file found for {run_accession}, skipping...")
            continue
        
        # Determine sample and PCR IDs
        bio_sample = sample_accession if sample_accession and sample_accession != '' else run_accession
        pcr_id = run_accession  # PCR ID is the run_accession
        
        # RGLB: use library_name if available, otherwise use sample_accession
        rglb = library_name if library_name and library_name != '' else sample_accession if sample_accession else run_accession
        
        samples_data.append({
            'run_accession': run_accession,
            'bio_sample': bio_sample,
            'pcr_id': pcr_id,
            'r1_url': r1_url,
            'r1_url_ftp': r1_url_ftp,
            'r2_url': r2_url,
            'r2_url_ftp': r2_url_ftp,
            'rglb': rglb,
            'library_name': library_name,
            'study_accession': study_accession,
            'library_layout': library_layout,
            'file_size_bytes': file_size
        })
    
    return samples_data, total_size_bytes

def download_fastq_files(samples_data, output_dir, dry_run=False, max_workers=8, use_aria2c=True, connections_per_file=2):
    """Download all FASTQ files in parallel."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    failed_runs = []
    has_aria2c = check_aria2c() if use_aria2c else False
    if use_aria2c and not has_aria2c:
        print("[WARNING] aria2c not found. Falling back to parallel wget/curl.")
        use_aria2c = False
    
    download_tasks = []
    sample_info_map = {}
    
    for sample_info in samples_data:
        run_accession = sample_info['run_accession']
        pcr_id = sample_info['pcr_id']
        r1_url = sample_info['r1_url']
        r2_url = sample_info.get('r2_url')
        
        r1_filename = os.path.basename(r1_url)
        r1_path = output_dir / r1_filename
        
        r2_path = None
        if r2_url:
            r2_filename = os.path.basename(r2_url)
            r2_path = output_dir / r2_filename
        
        sample_info['r1_path'] = str(r1_path)
        sample_info['r2_path'] = str(r2_path) if r2_path else ''
        
        # Use clearer task IDs: _SE for single-end, _R1/_R2 for paired-end
        library_layout = sample_info.get('library_layout', '').upper()
        if library_layout == 'SINGLE' or not r2_url:
            task_id = f"{run_accession}_SE"
        else:
            task_id = f"{run_accession}_R1"
        r1_ftp = sample_info.get('r1_url_ftp')
        download_tasks.append((task_id, r1_url, r1_path, sample_info, r1_ftp))
        sample_info_map[task_id] = sample_info
        
        if r2_url and r2_path:
            task_id = f"{run_accession}_R2"
            r2_ftp = sample_info.get('r2_url_ftp')
            download_tasks.append((task_id, r2_url, r2_path, sample_info, r2_ftp))
            sample_info_map[task_id] = sample_info
    
    if dry_run:
        print(f"\n[DRY RUN] Would download {len(download_tasks)} files...")
        print(f"Output directory: {output_dir}\n")
        for task_id, url, path, sample_info, ftp_url in download_tasks:
            print(f"  {task_id}: {url} -> {path}")
            if ftp_url:
                print(f"    (FTP fallback: {ftp_url})")
        return samples_data, []
    
    print(f"\nDownloading {len(download_tasks)} files using {'aria2c' if use_aria2c else f'{max_workers} parallel workers'}...")
    print(f"Output directory: {output_dir}\n")
    
    completed_samples = set()
    downloaded_files = []
    
    if use_aria2c:
        print(f"Using aria2c with {max_workers} parallel workers...\n")
        tasks_to_download = []
        for task_id, url, path, sample_info, ftp_url in download_tasks:
            if path.exists() and path.stat().st_size > 0:
                print(f"  ✓ {task_id}: Already exists")
                run_accession = sample_info['run_accession']
                if run_accession not in completed_samples:
                    r1_done = Path(sample_info['r1_path']).exists() and Path(sample_info['r1_path']).stat().st_size > 0
                    r2_done = not sample_info.get('r2_path') or (Path(sample_info['r2_path']).exists() and Path(sample_info['r2_path']).stat().st_size > 0)
                    if r1_done and r2_done:
                        completed_samples.add(run_accession)
                        downloaded_files.append(sample_info)
            else:
                tasks_to_download.append((task_id, url, path, sample_info, ftp_url))
        
        failed_tasks = []
        if tasks_to_download:
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                future_to_task = {}
                for task_id, url, path, sample_info, ftp_url in tasks_to_download:
                    future = executor.submit(
                        download_file_wrapper, 
                        url, path, use_aria2c=True, ftp_url=ftp_url, connections_per_file=connections_per_file
                    )
                    future_to_task[future] = (task_id, url, path, sample_info, ftp_url)
                
                completed_count = 0
                total_tasks = len(tasks_to_download)
                for future in as_completed(future_to_task):
                    task_id, url, path, sample_info, ftp_url = future_to_task[future]
                    completed_count += 1
                    try:
                        success, msg, error_msg = future.result()
                        if success:
                            print(f"  ✓ [{completed_count}/{total_tasks}] {task_id}: {msg}")
                        else:
                            print(f"  ✗ [{completed_count}/{total_tasks}] {task_id}: {msg}")
                            if error_msg:
                                print(f"      Error: {error_msg}")
                            failed_tasks.append((task_id, url, path, sample_info, ftp_url))
                            failed_runs.append((task_id, url, path, sample_info, ftp_url))
                    except Exception as e:
                        print(f"  ✗ [{completed_count}/{total_tasks}] {task_id}: Exception - {e}")
                        failed_tasks.append((task_id, url, path, sample_info, ftp_url))
                        failed_runs.append((task_id, url, path, sample_info, ftp_url))
                    
                    run_accession = sample_info['run_accession']
                    if run_accession not in completed_samples:
                        r1_done = Path(sample_info['r1_path']).exists() and Path(sample_info['r1_path']).stat().st_size > 0
                        r2_done = not sample_info.get('r2_path') or (Path(sample_info['r2_path']).exists() and Path(sample_info['r2_path']).stat().st_size > 0)
                        if r1_done and r2_done:
                            completed_samples.add(run_accession)
                            downloaded_files.append(sample_info)
    else:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future_to_task = {}
            for task_id, url, path, sample_info, ftp_url in download_tasks:
                if not path.exists() or path.stat().st_size == 0:
                    future = executor.submit(
                        download_file_wrapper, 
                        url, path, use_aria2c=False, ftp_url=ftp_url, connections_per_file=connections_per_file
                    )
                    future_to_task[future] = (task_id, url, path, sample_info, ftp_url)
            
            for future in as_completed(future_to_task):
                task_id, url, path, sample_info, ftp_url = future_to_task[future]
                try:
                    success, msg, error_msg = future.result()
                    if success:
                        print(f"  ✓ {task_id}: {msg}")
                    else:
                        print(f"  ✗ {task_id}: {msg}")
                        if error_msg:
                            print(f"      Error: {error_msg}")
                        failed_runs.append((task_id, url, path, sample_info, ftp_url))
                except Exception as e:
                    print(f"  ✗ {task_id}: Exception - {e}")
                    failed_runs.append((task_id, url, path, sample_info, ftp_url))
                
                run_accession = sample_info['run_accession']
                if run_accession not in completed_samples:
                    r1_done = Path(sample_info['r1_path']).exists() and Path(sample_info['r1_path']).stat().st_size > 0
                    r2_done = not sample_info.get('r2_path') or (Path(sample_info['r2_path']).exists() and Path(sample_info['r2_path']).stat().st_size > 0)
                    if r1_done and r2_done:
                        completed_samples.add(run_accession)
                        downloaded_files.append(sample_info)
    
    for sample_info in samples_data:
        run_accession = sample_info['run_accession']
        if run_accession not in completed_samples:
            r1_path = Path(sample_info['r1_path'])
            r2_path = Path(sample_info['r2_path']) if sample_info.get('r2_path') else None
            r1_done = r1_path.exists() and r1_path.stat().st_size > 0
            r2_done = not r2_path or (r2_path.exists() and r2_path.stat().st_size > 0)
            if r1_done and r2_done:
                completed_samples.add(run_accession)
                downloaded_files.append(sample_info)
    
    return downloaded_files, failed_runs

def create_samples_tsv(downloaded_files, output_tsv):
    """Create samples.tsv file in PIGSTI format."""
    rows = []
    
    for sample_info in downloaded_files:
        bio_sample = sample_info['bio_sample']
        pcr_id = sample_info['pcr_id']
        r1_path = sample_info['r1_path']
        r2_path = sample_info.get('r2_path', '')
        rglb = sample_info['rglb']
        study_accession = sample_info.get('study_accession', '')
        
        if not r2_path or r2_path == '':
            r2_path = ''
        
        rows.append({
            'sample': bio_sample,
            'pcr': pcr_id,
            'r1': r1_path,
            'r2': r2_path,
            'RGLB': rglb,
            'sequencing_run': study_accession if study_accession else 'ENA',
            # Mark ENA-derived rows (metadata only; pipeline does not delete raw FASTQs).
            'source': 'ENA'
        })
    
    df = pd.DataFrame(rows)
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"\n✓ Created samples.tsv: {output_tsv}")
    print(f"  Total samples: {len(df)}")
    print(f"  Unique biological samples: {df['sample'].nunique()}")
    print(f"  Total PCRs: {len(df)}")
    print(f"  Paired-end: {len(df[df['r2'] != ''])}")
    print(f"  Single-end: {len(df[df['r2'] == ''])}")
    
    return df

def main():
    parser = argparse.ArgumentParser(
        description='Filter ENA sheets for RANDOM library selection, download reads, and create PIGSTI samples.tsv',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('ena_files', nargs='+', help='ENA file report(s) (TSV format) - can specify multiple files')
    parser.add_argument('--output-dir', required=True, help='Output directory for FASTQ files')
    parser.add_argument('--samples-tsv', default=None,
                       help='Output samples.tsv file path (default: <output_dir>/samples.tsv)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be downloaded without actually downloading')
    parser.add_argument('--max-workers', type=int, default=8,
                       help='Maximum number of parallel downloads (default: 8)')
    parser.add_argument('--connections-per-file', type=int, default=2,
                       help='Number of connections per file for aria2c (default: 2)')
    parser.add_argument('--no-aria2c', action='store_true',
                       help='Disable aria2c and use wget/curl with multiprocessing instead')
    parser.add_argument('--no-filter-random', action='store_true',
                       help='Do not filter for library_selection=RANDOM (download all samples)')
    
    args = parser.parse_args()
    
    # Check if ENA files exist
    for ena_file in args.ena_files:
        if not os.path.exists(ena_file):
            print(f"ERROR: ENA file not found: {ena_file}")
            sys.exit(1)
    
    # Set default samples.tsv path
    if args.samples_tsv is None:
        args.samples_tsv = os.path.join(args.output_dir, 'samples.tsv')
    
    print("=" * 80)
    print("ENA to PIGSTI Sample Sheet Generator")
    print("=" * 80)
    print(f"ENA files: {', '.join(args.ena_files)}")
    print(f"Output directory: {args.output_dir}")
    print(f"Samples.tsv output: {args.samples_tsv}")
    print(f"Filter for RANDOM library selection: {not args.no_filter_random}")
    print("=" * 80)
    
    # Parse all ENA files
    print("\n[1/3] Parsing ENA files...")
    all_samples_data = []
    total_size_all_files = 0
    samples_with_size = 0
    samples_without_size = 0
    
    for ena_file in args.ena_files:
        print(f"\nProcessing: {ena_file}")
        samples_data, total_size = parse_ena_file(ena_file, filter_random=not args.no_filter_random)
        print(f"  Found {len(samples_data)} samples with valid FTP URLs")
        all_samples_data.extend(samples_data)
        total_size_all_files += total_size
        
        # Count samples with/without size info
        for sample in samples_data:
            if sample.get('file_size_bytes', 0) > 0:
                samples_with_size += 1
            else:
                samples_without_size += 1
    
    if len(all_samples_data) == 0:
        print("ERROR: No valid samples found in ENA files!")
        sys.exit(1)
    
    print(f"\nTotal samples across all files: {len(all_samples_data)}")
    
    # Calculate download statistics
    paired_count = sum(1 for s in all_samples_data if s.get('r2_url'))
    single_count = len(all_samples_data) - paired_count
    total_files = len(all_samples_data) + paired_count  # R1 for all + R2 for paired
    
    # Display size estimate and ask for confirmation
    print("\n" + "=" * 80)
    print("DOWNLOAD SIZE ESTIMATE")
    print("=" * 80)
    print(f"Total samples: {len(all_samples_data)}")
    print(f"  - Paired-end: {paired_count}")
    print(f"  - Single-end: {single_count}")
    print(f"Total files to download: {total_files}")
    print(f"\nEstimated total size: {format_bytes(total_size_all_files)}")
    
    if samples_without_size > 0:
        print(f"\n[NOTE] Size information available for {samples_with_size} samples")
        print(f"        Size information missing for {samples_without_size} samples")
        print(f"        The actual download size may be larger than estimated.")
    
    print(f"\nOutput directory: {args.output_dir}")
    print("=" * 80)
    
    # Ask for user confirmation
    if not args.dry_run:
        while True:
            response = input("\nProceed with download? (yes/no): ").strip().lower()
            if response in ['yes', 'y']:
                print("\nStarting downloads...")
                break
            elif response in ['no', 'n']:
                print("\nDownload cancelled by user.")
                sys.exit(0)
            else:
                print("Please enter 'yes' or 'no'")
    else:
        print("\n[DRY RUN] Skipping download confirmation (dry-run mode)")
    
    # Download files
    print("\n[2/3] Downloading FASTQ files...")
    downloaded_files, failed_runs = download_fastq_files(
        all_samples_data, 
        args.output_dir, 
        dry_run=args.dry_run,
        max_workers=args.max_workers,
        use_aria2c=not args.no_aria2c,
        connections_per_file=args.connections_per_file
    )
    
    if len(downloaded_files) == 0 and not args.dry_run:
        print("ERROR: No files were downloaded!")
        sys.exit(1)
    
    # Retry failed downloads with reduced parallelism and delays (to avoid rate limiting)
    if failed_runs and not args.dry_run:
        # Track which samples are already complete
        completed_samples = set()
        for sample_info in downloaded_files:
            completed_samples.add(sample_info['run_accession'])
        print(f"\n[RETRY] Retrying {len(failed_runs)} failed downloads with reduced parallelism...")
        print("  (Using fewer workers and delays to avoid ENA server rate limiting)")
        print("")
        
        # Use fewer workers for retry (max 4, or 1/4 of original)
        retry_workers = min(4, max(1, args.max_workers // 4))
        retry_tasks = failed_runs.copy()
        failed_runs = []  # Clear for retry round
        
        # Add small delay before retry to let server recover
        time.sleep(5)
        
        has_aria2c = check_aria2c() if not args.no_aria2c else False
        use_aria2c = not args.no_aria2c and has_aria2c
        
        if use_aria2c:
            # Retry with aria2c but fewer connections and workers
            with ThreadPoolExecutor(max_workers=retry_workers) as executor:
                future_to_task = {}
                for i, (task_id, url, path, sample_info, ftp_url) in enumerate(retry_tasks):
                    # Stagger the retries with small delays
                    if i > 0:
                        time.sleep(2)
                    
                    future = executor.submit(
                        download_file_wrapper, 
                        url, path, use_aria2c=True, ftp_url=ftp_url, connections_per_file=1  # Use only 1 connection per file on retry
                    )
                    future_to_task[future] = (task_id, url, path, sample_info, ftp_url)
                
                completed_count = 0
                total_retries = len(retry_tasks)
                for future in as_completed(future_to_task):
                    task_id, url, path, sample_info, ftp_url = future_to_task[future]
                    completed_count += 1
                    try:
                        success, msg, error_msg = future.result()
                        if success:
                            print(f"  ✓ [RETRY {completed_count}/{total_retries}] {task_id}: {msg}")
                            # Update downloaded_files if this completes a sample
                            run_accession = sample_info['run_accession']
                            if run_accession not in completed_samples:
                                r1_done = Path(sample_info['r1_path']).exists() and Path(sample_info['r1_path']).stat().st_size > 0
                                r2_done = not sample_info.get('r2_path') or (Path(sample_info['r2_path']).exists() and Path(sample_info['r2_path']).stat().st_size > 0)
                                if r1_done and r2_done:
                                    completed_samples.add(run_accession)
                                    downloaded_files.append(sample_info)
                        else:
                            print(f"  ✗ [RETRY {completed_count}/{total_retries}] {task_id}: {msg}")
                            if error_msg:
                                print(f"      Error: {error_msg[:100]}")
                            failed_runs.append((task_id, url, path, sample_info, ftp_url))
                    except Exception as e:
                        print(f"  ✗ [RETRY {completed_count}/{total_retries}] {task_id}: Exception - {e}")
                        failed_runs.append((task_id, url, path, sample_info, ftp_url))
        else:
            # Retry with wget/curl, single-threaded with delays
            for i, (task_id, url, path, sample_info, ftp_url) in enumerate(retry_tasks):
                if i > 0:
                    time.sleep(3)  # Delay between retries
                
                success, msg, error_msg = download_file_wrapper(
                    url, path, use_aria2c=False, ftp_url=ftp_url, connections_per_file=1
                )
                if success:
                    print(f"  ✓ [RETRY {i+1}/{len(retry_tasks)}] {task_id}: {msg}")
                    run_accession = sample_info['run_accession']
                    if run_accession not in completed_samples:
                        r1_done = Path(sample_info['r1_path']).exists() and Path(sample_info['r1_path']).stat().st_size > 0
                        r2_done = not sample_info.get('r2_path') or (Path(sample_info['r2_path']).exists() and Path(sample_info['r2_path']).stat().st_size > 0)
                        if r1_done and r2_done:
                            completed_samples.add(run_accession)
                            downloaded_files.append(sample_info)
                else:
                    print(f"  ✗ [RETRY {i+1}/{len(retry_tasks)}] {task_id}: {msg}")
                    if error_msg:
                        print(f"      Error: {error_msg[:100]}")
                    failed_runs.append((task_id, url, path, sample_info, ftp_url))
        
        if failed_runs:
            print(f"\n⚠ {len(failed_runs)} files still failed after retry")
        else:
            print(f"\n✓ All files successfully downloaded after retry!")
    
    # Write failed runs list if any remain
    if failed_runs:
        failed_file = Path(args.output_dir) / 'failed_downloads.tsv'
        rows = []
        for task_id, url, path, sample_info, ftp_url in failed_runs:
            rows.append({
                'run_accession': sample_info['run_accession'],
                'bio_sample': sample_info.get('bio_sample', ''),
                'pcr_id': sample_info.get('pcr_id', ''),
                'file_type': 'R1' if '_R1' in task_id else 'R2' if '_R2' in task_id else 'SE',
                'http_url': url,
                'ftp_url': ftp_url if ftp_url else '',
                'expected_filename': Path(path).name,
                'expected_path': str(path)
            })
        df_failed = pd.DataFrame(rows)
        df_failed.to_csv(failed_file, sep='\t', index=False)
        print(f"\n⚠ Created failed downloads list: {failed_file}")
        print(f"  Total failed files: {len(df_failed)}")
        print(f"  You can manually download these files or re-run the script later.")
    
    # Create samples.tsv
    print("\n[3/3] Creating samples.tsv...")
    df = create_samples_tsv(downloaded_files, args.samples_tsv)
    
    print("\n" + "=" * 80)
    print("✓ Complete!")
    print("=" * 80)
    print(f"\nNext steps:")
    print(f"  1. Review the samples.tsv file: {args.samples_tsv}")
    if failed_runs:
        print(f"  2. Check failed_downloads.tsv in {args.output_dir} for failed downloads")
        print(f"  3. You can manually download failed files or re-run the script")
    else:
        print(f"  2. Copy it to your PIGSTI config directory if needed")
        print(f"  3. Update the config.yaml to point to the correct samples.tsv path")
    print()

if __name__ == '__main__':
    main()

