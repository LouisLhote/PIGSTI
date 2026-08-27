#!/usr/bin/env python3
"""
Download FASTQ files from ENA and create samples.tsv file.

Usage:
    python download_ena_fastq.py <ena_file> <output_dir> [--samples-tsv <output_tsv>]

Example:
    python download_ena_fastq.py filereport_read_run_PRJEB105076.tsv /raid_md0/louis/reanalysis/Miranda_2026/raw_datas/
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
from functools import partial

def check_aria2c():
    """Check if aria2c is available."""
    try:
        subprocess.run(['aria2c', '--version'], capture_output=True, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False

def download_file_aria2c(url, output_path, max_connections=4, max_retries=5, retry_wait=3):
    """Download a file using aria2c (parallel downloader)."""
    # aria2c options optimized for ENA FTP servers:
    # -x: max connections per server
    # -s: split into N connections
    # --max-tries: max retry attempts
    # --retry-wait: wait time between retries
    # --timeout: connection timeout (longer for large files)
    # --connect-timeout: connection establishment timeout
    # --min-split-size: minimum size to split (avoid splitting small files)
    cmd = [
        'aria2c',
        '-x', str(max_connections),  # Max connections per server
        '-s', str(max_connections),  # Split into N connections
        '--min-split-size', '10M',  # Only split files > 10MB
        '--max-tries', str(max_retries),
        '--retry-wait', str(retry_wait),
        '--timeout', '60',  # Connection timeout (increased for large files)
        '--connect-timeout', '20',  # Connection establishment timeout
        '--continue',  # Resume partial downloads
        '--file-allocation', 'none',  # Don't pre-allocate (faster for large files)
        '--allow-overwrite',  # Overwrite existing files if needed
        '--auto-file-renaming=false',  # Don't rename on conflict
        '--max-connection-per-server', str(max_connections),  # Ensure this matches -x
        '--split', str(max_connections),  # Ensure this matches -s
        '--console-log-level=warn',  # Show warnings and errors (not quiet, so we can see errors)
        '--summary-interval=0',  # Don't show progress summary
        '-d', str(output_path.parent),
        '-o', output_path.name,
        url
    ]
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=10800)  # 3 hour timeout
        # Verify file was actually downloaded
        if output_path.exists() and output_path.stat().st_size > 0:
            return True, None
        else:
            return False, "File downloaded but is empty or missing"
    except subprocess.TimeoutExpired:
        return False, "Download timeout (exceeded 3 hours)"
    except subprocess.CalledProcessError as e:
        # Capture both stderr and stdout for better error reporting
        error_parts = []
        if e.stderr and e.stderr.strip():
            error_parts.append(f"stderr: {e.stderr.strip()}")
        if e.stdout and e.stdout.strip():
            error_parts.append(f"stdout: {e.stdout.strip()}")
        if not error_parts:
            # Try to get error from aria2c log or check if file exists
            if not output_path.exists():
                error_parts.append(f"File not created (exit code {e.returncode})")
            else:
                error_parts.append(f"Exit code: {e.returncode}")
        
        error_msg = " | ".join(error_parts) if error_parts else f"Exit code: {e.returncode}"
        
        # Extract meaningful error message
        error_lower = error_msg.lower()
        if "not found" in error_lower or "404" in error_msg or "HTTP/1.1 404" in error_msg:
            return False, "File not found (404)"
        elif "timeout" in error_lower or "timed out" in error_lower:
            return False, "Connection timeout"
        elif "connection" in error_lower or "refused" in error_lower or "could not resolve" in error_lower:
            return False, "Connection error"
        elif "403" in error_msg or "Forbidden" in error_msg:
            return False, "Access forbidden (403)"
        elif "503" in error_msg or "service unavailable" in error_lower:
            return False, "Service unavailable (503)"
        elif e.returncode == 1:
            return False, f"aria2c failed (exit code 1) - check network/server"
        else:
            # Return more detailed error
            return False, f"aria2c failed (code {e.returncode}): {error_msg[:200]}"

def download_file_single(url, output_path, max_retries=5, verbose=False, ftp_url=None):
    """Download a file from URL with retries (single-threaded)."""
    last_error = None
    urls_to_try = [url]
    if ftp_url:
        urls_to_try.append(ftp_url)
    
    for attempt in range(max_retries):
        # Progressive backoff: longer delays on later attempts
        if attempt > 0:
            delay = min(10 * attempt, 30)  # 10s, 20s, 30s, 30s, 30s
            if verbose:
                print(f"        Retry {attempt + 1}/{max_retries} after {delay}s delay...")
            time.sleep(delay)
        
        # Try each URL (HTTP first, then FTP if available)
        for current_url in urls_to_try:
            try:
                # Try wget first
                # For ENA, ftp.sra.ebi.ac.uk URLs work better with HTTP protocol
                # When you run wget manually without protocol, it defaults to HTTP
                if current_url.startswith('ftp://'):
                    # For FTP URLs, try HTTP first (ENA servers support both)
                    # Convert ftp:// to http:// for better reliability (matches manual wget behavior)
                    http_url = current_url.replace('ftp://', 'http://', 1)
                    cmd = ['wget', '-O', str(output_path), '--continue', '--timeout=120', '--tries=1', '--waitretry=5', http_url]
                elif current_url.startswith('http://') or current_url.startswith('https://'):
                    # Already HTTP/HTTPS
                    cmd = ['wget', '-O', str(output_path), '--continue', '--timeout=120', '--tries=1', '--waitretry=5', current_url]
                else:
                    # No protocol specified - explicitly use HTTP (matches manual wget behavior)
                    http_url = 'http://' + current_url
                    cmd = ['wget', '-O', str(output_path), '--continue', '--timeout=120', '--tries=1', '--waitretry=5', http_url]
                
                result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=10800)  # 3 hour timeout
                # Verify file was downloaded
                if output_path.exists() and output_path.stat().st_size > 0:
                    return True, None
                else:
                    # File exists but is empty or missing - try next URL or retry
                    if output_path.exists():
                        output_path.unlink()  # Remove empty file
                    last_error = "File downloaded but is empty"
                    continue
            except subprocess.TimeoutExpired:
                # Clean up partial file on timeout
                if output_path.exists():
                    try:
                        output_path.unlink()
                    except:
                        pass
                last_error = "Download timeout"
                continue
            except subprocess.CalledProcessError as e:
                # Clean up partial file on error
                if output_path.exists():
                    try:
                        output_path.unlink()
                    except:
                        pass
                
                # Capture error message
                error_msg = ""
                if e.stderr:
                    error_msg = e.stderr.strip()
                elif e.stdout:
                    error_msg = e.stdout.strip()
                
                # Check for specific error codes
                if "404" in error_msg or "Not Found" in error_msg:
                    # If this was HTTP and we have FTP, try FTP next
                    if current_url.startswith('http') and ftp_url and current_url == urls_to_try[0]:
                        continue  # Try FTP URL
                    return False, "File not found (404)"
                elif "403" in error_msg or "Forbidden" in error_msg:
                    return False, "Access forbidden (403)"
                elif "Connection refused" in error_msg or "could not resolve" in error_msg.lower() or "Connection reset" in error_msg or "failed: Connection refused" in error_msg:
                    # If this was HTTP and we have FTP, try FTP next
                    if current_url.startswith('http') and ftp_url and current_url == urls_to_try[0]:
                        last_error = "HTTP connection failed, will try FTP"
                        continue  # Try FTP URL
                    # Also try FTP if we haven't tried it yet
                    if ftp_url and current_url != ftp_url and ftp_url not in [u for u in urls_to_try if u != current_url]:
                        last_error = "HTTP connection failed, will try FTP"
                        continue  # Try FTP URL
                    last_error = f"Connection error: {error_msg[:100]}" if error_msg else "Connection error"
                    continue
                else:
                    # If this was HTTP and we have FTP, try FTP next
                    if current_url.startswith('http') and ftp_url and current_url == urls_to_try[0]:
                        last_error = f"HTTP failed: {error_msg[:100] if error_msg else 'Unknown error'}, will try FTP"
                        continue  # Try FTP URL
                    last_error = f"wget error: {error_msg[:100]}" if error_msg else f"wget failed (exit code {e.returncode})"
                    continue
            except FileNotFoundError:
                # Try curl if wget is not available
                for curl_url in urls_to_try:
                    try:
                        if curl_url.startswith('ftp://'):
                            cmd = ['curl', '-o', str(output_path), '-C', '-', '--max-time', '10800', '--retry', '3', '--retry-delay', '5', '--ftp-pasv', curl_url]
                        else:
                            cmd = ['curl', '-L', '-o', str(output_path), '-C', '-', '--max-time', '10800', '--retry', '3', '--retry-delay', '5', '--fail', curl_url]
                        result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=10800)
                        if output_path.exists() and output_path.stat().st_size > 0:
                            return True, None
                        last_error = "curl: File downloaded but is empty"
                    except subprocess.CalledProcessError as e:
                        error_msg = e.stderr.strip() if e.stderr else e.stdout.strip() if e.stdout else ""
                        if "404" in error_msg or "Not Found" in error_msg:
                            if curl_url == urls_to_try[0] and len(urls_to_try) > 1:
                                continue  # Try next URL
                            return False, "File not found (404)"
                        last_error = f"curl error: {error_msg[:100]}" if error_msg else f"curl failed (exit code {e.returncode})"
                        if curl_url == urls_to_try[0] and len(urls_to_try) > 1:
                            continue  # Try next URL
                    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
                        last_error = "Neither wget nor curl available" if isinstance(e, FileNotFoundError) else "Download timeout"
                        break
                continue
    
    return False, last_error or "Download failed after all retries"

def download_file_wrapper(url, output_path, use_aria2c=False, max_workers=1, retry_count=0, connections_per_file=2, ftp_url=None):
    """Wrapper function for downloading a single file."""
    if output_path.exists() and output_path.stat().st_size > 0:
        return (True, f"Already exists: {output_path.name}", None)
    
    if use_aria2c:
        # Reduce connections on retry
        conns = max(1, connections_per_file - retry_count)
        success, error_msg = download_file_aria2c(url, output_path, max_connections=conns, max_retries=5, retry_wait=5)
        
        # If aria2c fails, fallback to wget/curl (only on first attempt, not retries)
        if not success and retry_count == 0:
            # Add a small delay before fallback to avoid hammering the server
            time.sleep(2)
            # Try once with wget/curl as fallback (with more retries, and FTP if available)
            print(f"      Trying fallback (wget/curl) for {output_path.name}...")
            fallback_success, fallback_error = download_file_single(url, output_path, max_retries=5, verbose=True, ftp_url=ftp_url)
            if fallback_success:
                return (True, f"Downloaded (fallback): {output_path.name}", None)
            else:
                combined_error = f"aria2c: {error_msg}"
                if fallback_error:
                    combined_error += f"; wget/curl: {fallback_error}"
                return (False, f"Failed: {output_path.name}", combined_error)
    else:
        success, error_msg = download_file_single(url, output_path, max_retries=5, ftp_url=ftp_url)
        if not success:
            error_msg = error_msg or "Download failed"
    
    if success:
        return (True, f"Downloaded: {output_path.name}", None)
    else:
        return (False, f"Failed: {output_path.name}", error_msg)

def parse_ena_file(ena_file):
    """Parse ENA file report and extract download information."""
    df = pd.read_csv(ena_file, sep='\t')
    
    samples_data = []
    
    for idx, row in df.iterrows():
        run_accession = row['run_accession']
        sample_accession = row.get('sample_accession', '')
        study_accession = row.get('study_accession', '')
        library_name = row.get('library_name', '')
        submitted_ftp = row.get('submitted_ftp', '')
        
        # Skip if no FTP URLs
        if pd.isna(submitted_ftp) or submitted_ftp == '':
            print(f"[WARNING] No FTP URLs for {run_accession}, skipping...")
            continue
        
        # Parse FTP URLs (semicolon-separated)
        ftp_urls = [url.strip() for url in str(submitted_ftp).split(';') if url.strip()]
        
        # Find R1 and R2 files
        r1_url = None
        r1_url_ftp = None
        r2_url = None
        r2_url_ftp = None
        
        for url in ftp_urls:
            # Keep original URL
            original_url = url.strip()
            
            # Create HTTP version (try HTTP first, fallback to FTP)
            if original_url.startswith('ftp://'):
                http_url = original_url.replace('ftp://', 'http://', 1)
                ftp_url = original_url
            elif not original_url.startswith('http'):
                # No protocol specified, assume FTP
                http_url = 'http://' + original_url
                ftp_url = 'ftp://' + original_url
            else:
                # Already HTTP
                http_url = original_url
                ftp_url = original_url.replace('http://', 'ftp://', 1) if original_url.startswith('http://') else None
            
            filename = os.path.basename(original_url)
            # Store both HTTP and FTP URLs
            if '_R1_' in filename or '_r1_' in filename:
                r1_url = http_url
                r1_url_ftp = ftp_url
            elif '_R2_' in filename or '_r2_' in filename:
                r2_url = http_url
                r2_url_ftp = ftp_url
        
        if not r1_url:
            print(f"[WARNING] No R1 file found for {run_accession}, skipping...")
            continue
        
        # Determine sample and PCR IDs based on user's correction:
        # bio_sample = sample_accession
        # pcr_id = run_accession
        # sequencing_run = study_accession
        
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
            'study_accession': study_accession
        })
    
    return samples_data

def download_fastq_files(samples_data, output_dir, dry_run=False, max_workers=8, use_aria2c=True, connections_per_file=2):
    """Download all FASTQ files in parallel."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Track failed downloads - initialize early
    failed_runs = []
    
    # Check for aria2c
    has_aria2c = check_aria2c() if use_aria2c else False
    if use_aria2c and not has_aria2c:
        print("[WARNING] aria2c not found. Falling back to parallel wget/curl.")
        use_aria2c = False
    
    # Prepare download tasks
    download_tasks = []
    sample_info_map = {}
    
    for sample_info in samples_data:
        run_accession = sample_info['run_accession']
        pcr_id = sample_info['pcr_id']
        r1_url = sample_info['r1_url']
        r2_url = sample_info.get('r2_url')
        
        # Determine output filenames
        r1_filename = os.path.basename(r1_url)
        r1_path = output_dir / r1_filename
        
        r2_path = None
        if r2_url:
            r2_filename = os.path.basename(r2_url)
            r2_path = output_dir / r2_filename
        
        # Store paths
        sample_info['r1_path'] = str(r1_path)
        sample_info['r2_path'] = str(r2_path) if r2_path else ''
        
        # Add download tasks (with FTP fallback URLs)
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
        for task_item in download_tasks:
            if len(task_item) == 5:
                task_id, url, path, sample_info, ftp_url = task_item
                print(f"  {task_id}: {url} -> {path}")
                if ftp_url:
                    print(f"    (FTP fallback: {ftp_url})")
            else:
                task_id, url, path, sample_info = task_item[:4]
                print(f"  {task_id}: {url} -> {path}")
        return samples_data
    
    print(f"\nDownloading {len(download_tasks)} files using {'aria2c' if use_aria2c else f'{max_workers} parallel workers'}...")
    print(f"Output directory: {output_dir}\n")
    
    # Track completed samples
    completed_samples = set()
    downloaded_files = []
    
    if use_aria2c:
        # Use aria2c with parallel processes (multiple aria2c instances running simultaneously)
        print(f"Using aria2c with {max_workers} parallel workers...\n")
        
        # Filter out files that already exist
        tasks_to_download = []
        for task_item in download_tasks:
            if len(task_item) == 5:
                task_id, url, path, sample_info, ftp_url = task_item
            else:
                task_id, url, path, sample_info = task_item[:4]
                ftp_url = None
            
            if path.exists() and path.stat().st_size > 0:
                print(f"  ✓ {task_id}: Already exists")
                # Still track completion
                run_accession = sample_info['run_accession']
                if run_accession not in completed_samples:
                    r1_done = Path(sample_info['r1_path']).exists() and Path(sample_info['r1_path']).stat().st_size > 0
                    r2_done = not sample_info.get('r2_path') or (Path(sample_info['r2_path']).exists() and Path(sample_info['r2_path']).stat().st_size > 0)
                    if r1_done and r2_done:
                        completed_samples.add(run_accession)
                        downloaded_files.append(sample_info)
            else:
                tasks_to_download.append((task_id, url, path, sample_info, ftp_url))
        
        # Download remaining files in parallel using ThreadPoolExecutor
        failed_tasks = []
        if tasks_to_download:
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                # Submit all download tasks
                future_to_task = {}
                for task_item in tasks_to_download:
                    if len(task_item) == 5:
                        task_id, url, path, sample_info, ftp_url = task_item
                    else:
                        task_id, url, path, sample_info = task_item[:4]
                        ftp_url = None
                    
                    future = executor.submit(
                        download_file_wrapper, 
                        url, path, use_aria2c=True, max_workers=1, retry_count=0, connections_per_file=connections_per_file, ftp_url=ftp_url
                    )
                    future_to_task[future] = (task_id, url, path, sample_info, ftp_url)
                
                # Process completed downloads with progress tracking
                completed_count = 0
                total_tasks = len(tasks_to_download)
                for future in as_completed(future_to_task):
                    task_item = future_to_task[future]
                    if len(task_item) == 5:
                        task_id, url, path, sample_info, ftp_url = task_item
                    else:
                        task_id, url, path, sample_info = task_item[:4]
                        ftp_url = None
                    
                    completed_count += 1
                    try:
                        success, msg, error_msg = future.result()
                        if success:
                            print(f"  ✓ [{completed_count}/{total_tasks}] {task_id}: {msg}")
                        else:
                            print(f"  ✗ [{completed_count}/{total_tasks}] {task_id}: {msg}")
                            if error_msg:
                                print(f"      Error: {error_msg}")
                            # Store for retry
                            failed_item = (task_id, url, path, sample_info, ftp_url)
                            failed_tasks.append(failed_item)
                            failed_runs.append(failed_item)
                    except Exception as e:
                        print(f"  ✗ [{completed_count}/{total_tasks}] {task_id}: Exception - {e}")
                        failed_item = (task_id, url, path, sample_info, ftp_url)
                        failed_tasks.append(failed_item)
                        failed_runs.append(failed_item)
                    
                    # Track completed samples
                    run_accession = sample_info['run_accession']
                    if run_accession not in completed_samples:
                        # Check if both R1 and R2 are done (if R2 exists)
                        r1_done = Path(sample_info['r1_path']).exists() and Path(sample_info['r1_path']).stat().st_size > 0
                        r2_done = not sample_info.get('r2_path') or (Path(sample_info['r2_path']).exists() and Path(sample_info['r2_path']).stat().st_size > 0)
                        if r1_done and r2_done:
                            completed_samples.add(run_accession)
                            downloaded_files.append(sample_info)
        
        # Retry failed downloads with reduced parallelism and longer waits
        if failed_tasks:
            print(f"\n[RETRY] Retrying {len(failed_tasks)} failed downloads with reduced parallelism...")
            print("  (Using wget/curl only for better reliability)\n")
            
            # Use fewer workers for retries to avoid overwhelming the server
            retry_workers = min(4, max_workers // 2)
            retry_tasks = failed_tasks.copy()
            failed_tasks = []
            
            with ThreadPoolExecutor(max_workers=retry_workers) as executor:
                future_to_task = {}
                for i, task_item in enumerate(retry_tasks):
                    if len(task_item) == 5:
                        task_id, url, path, sample_info, ftp_url = task_item
                    else:
                        task_id, url, path, sample_info = task_item[:4]
                        ftp_url = None
                    
                    # Add progressive delay before retrying to avoid rate limiting
                    time.sleep(i * 2)  # 0s, 2s, 4s, 6s, etc.
                    # Use wget/curl only for retries (more reliable)
                    future = executor.submit(
                        download_file_wrapper, 
                        url, path, use_aria2c=False, max_workers=1, retry_count=1, connections_per_file=connections_per_file, ftp_url=ftp_url
                    )
                    future_to_task[future] = (task_id, url, path, sample_info, ftp_url)
                
                completed_count = 0
                total_retries = len(retry_tasks)
                for future in as_completed(future_to_task):
                    task_item = future_to_task[future]
                    if len(task_item) == 5:
                        task_id, url, path, sample_info, ftp_url = task_item
                    else:
                        task_id, url, path, sample_info = task_item[:4]
                        ftp_url = None
                    
                    completed_count += 1
                    try:
                        success, msg, error_msg = future.result()
                        if success:
                            print(f"  ✓ [RETRY {completed_count}/{total_retries}] {task_id}: {msg}")
                        else:
                            print(f"  ✗ [RETRY {completed_count}/{total_retries}] {task_id}: {msg}")
                            if error_msg:
                                print(f"      Error: {error_msg}")
                            failed_item = (task_id, url, path, sample_info, ftp_url)
                            failed_tasks.append(failed_item)
                            # Also add to failed_runs if not already there
                            if failed_item not in failed_runs:
                                failed_runs.append(failed_item)
                    except Exception as e:
                        print(f"  ✗ [RETRY {completed_count}/{total_retries}] {task_id}: Exception - {e}")
                        failed_item = (task_id, url, path, sample_info, ftp_url)
                        failed_tasks.append(failed_item)
                        # Also add to failed_runs if not already there
                        if failed_item not in failed_runs:
                            failed_runs.append(failed_item)
                    
                    # Track completed samples
                    run_accession = sample_info['run_accession']
                    if run_accession not in completed_samples:
                        r1_done = Path(sample_info['r1_path']).exists() and Path(sample_info['r1_path']).stat().st_size > 0
                        r2_done = not sample_info.get('r2_path') or (Path(sample_info['r2_path']).exists() and Path(sample_info['r2_path']).stat().st_size > 0)
                        if r1_done and r2_done:
                            completed_samples.add(run_accession)
                            downloaded_files.append(sample_info)
            
            if failed_tasks:
                print(f"\n[WARNING] {len(failed_tasks)} files still failed after retry:")
                for task_item in failed_tasks:
                    if len(task_item) == 5:
                        task_id, url, path, sample_info, ftp_url = task_item
                    else:
                        task_id, url, path, sample_info = task_item[:4]
                        ftp_url = None
                    print(f"  - {task_id}: {url}")
                    if ftp_url:
                        print(f"    FTP alternative: {ftp_url}")
                    # Add to failed_runs list if not already there
                    if task_item not in failed_runs:
                        failed_runs.append(task_item)
                print("\nYou may need to download these manually or check the URLs.")
        
        # Add any remaining samples that already had all files
        for sample_info in samples_data:
            run_accession = sample_info['run_accession']
            if run_accession not in completed_samples:
                r1_done = Path(sample_info['r1_path']).exists()
                r2_done = not sample_info.get('r2_path') or Path(sample_info['r2_path']).exists()
                if r1_done and r2_done:
                    completed_samples.add(run_accession)
                    downloaded_files.append(sample_info)
    else:
        # Use ThreadPoolExecutor for parallel downloads (fallback when aria2c not available)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all download tasks
            future_to_task = {}
            for task_item in download_tasks:
                if len(task_item) == 5:
                    task_id, url, path, sample_info, ftp_url = task_item
                else:
                    task_id, url, path, sample_info = task_item[:4]
                    ftp_url = None
                
                if not path.exists() or path.stat().st_size == 0:
                    future = executor.submit(
                        download_file_wrapper, 
                        url, path, use_aria2c=False, max_workers=1, ftp_url=ftp_url
                    )
                    future_to_task[future] = (task_id, url, path, sample_info, ftp_url)
            
            # Process completed downloads
            for future in as_completed(future_to_task):
                task_item = future_to_task[future]
                if len(task_item) == 5:
                    task_id, url, path, sample_info, ftp_url = task_item
                else:
                    task_id, url, path, sample_info = task_item[:4]
                    ftp_url = None
                
                try:
                    success, msg, error_msg = future.result()
                    if success:
                        print(f"  ✓ {task_id}: {msg}")
                    else:
                        print(f"  ✗ {task_id}: {msg}")
                        if error_msg:
                            print(f"      Error: {error_msg}")
                        # Track failed downloads
                        failed_item = (task_id, url, path, sample_info, ftp_url)
                        if failed_item not in failed_runs:
                            failed_runs.append(failed_item)
                except Exception as e:
                    print(f"  ✗ {task_id}: Exception - {e}")
                    # Track failed downloads
                    failed_item = (task_id, url, path, sample_info, ftp_url)
                    if failed_item not in failed_runs:
                        failed_runs.append(failed_item)
                
                # Track completed samples
                run_accession = sample_info['run_accession']
                if run_accession not in completed_samples:
                    # Check if both R1 and R2 are done (if R2 exists)
                    r1_done = Path(sample_info['r1_path']).exists() and Path(sample_info['r1_path']).stat().st_size > 0
                    r2_done = not sample_info.get('r2_path') or (Path(sample_info['r2_path']).exists() and Path(sample_info['r2_path']).stat().st_size > 0)
                    if r1_done and r2_done:
                        completed_samples.add(run_accession)
                        downloaded_files.append(sample_info)
        
        # Add samples that already had all files
        for sample_info in samples_data:
            run_accession = sample_info['run_accession']
            if run_accession not in completed_samples:
                r1_done = Path(sample_info['r1_path']).exists() and Path(sample_info['r1_path']).stat().st_size > 0
                r2_done = not sample_info.get('r2_path') or (Path(sample_info['r2_path']).exists() and Path(sample_info['r2_path']).stat().st_size > 0)
                if r1_done and r2_done:
                    completed_samples.add(run_accession)
                    downloaded_files.append(sample_info)
    
    # Check for any samples that didn't complete (missing files)
    for sample_info in samples_data:
        run_accession = sample_info['run_accession']
        if run_accession not in completed_samples:
            r1_path = Path(sample_info['r1_path'])
            r2_path = Path(sample_info['r2_path']) if sample_info.get('r2_path') else None
            
            # Check if files are missing or empty
            r1_missing = not r1_path.exists() or r1_path.stat().st_size == 0
            r2_missing = r2_path and (not r2_path.exists() or r2_path.stat().st_size == 0)
            
            if r1_missing:
                task_id = f"{run_accession}_R1"
                url = sample_info['r1_url']
                ftp_url = sample_info.get('r1_url_ftp')
                failed_item = (task_id, url, r1_path, sample_info, ftp_url)
                if failed_item not in failed_runs:
                    failed_runs.append(failed_item)
            
            if r2_missing and r2_path:
                task_id = f"{run_accession}_R2"
                url = sample_info.get('r2_url', '')
                ftp_url = sample_info.get('r2_url_ftp')
                failed_item = (task_id, url, r2_path, sample_info, ftp_url)
                if failed_item not in failed_runs:
                    failed_runs.append(failed_item)
    
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
        
        # If single-end (no R2), use empty string
        if not r2_path or r2_path == '':
            r2_path = ''
        
        rows.append({
            'sample': bio_sample,
            'pcr': pcr_id,
            'r1': r1_path,
            'r2': r2_path,
            'RGLB': rglb,
            'sequencing_run': sample_info.get('study_accession', 'ENA_PRJEB105076')  # Study accession from ENA file
        })
    
    # Create DataFrame and save
    df = pd.DataFrame(rows)
    df.to_csv(output_tsv, sep='\t', index=False)
    print(f"\n✓ Created samples.tsv: {output_tsv}")
    print(f"  Total samples: {len(df)}")
    print(f"  Unique biological samples: {df['sample'].nunique()}")
    print(f"  Total PCRs: {len(df)}")
    
    return df

def write_failed_runs(failed_runs, output_dir):
    """Write failed runs to a TSV file for easy review and re-download."""
    if not failed_runs:
        return None
    
    failed_file = Path(output_dir) / 'failed_downloads.tsv'
    
    rows = []
    for failed_item in failed_runs:
        if len(failed_item) == 5:
            task_id, url, path, sample_info, ftp_url = failed_item
        else:
            task_id, url, path, sample_info = failed_item[:4]
            ftp_url = None
        
        run_accession = sample_info['run_accession']
        bio_sample = sample_info.get('bio_sample', run_accession)
        pcr_id = sample_info.get('pcr_id', run_accession)
        
        # Determine if R1 or R2
        file_type = 'R1' if '_R1' in task_id or 'R1' in task_id else 'R2'
        
        rows.append({
            'run_accession': run_accession,
            'bio_sample': bio_sample,
            'pcr_id': pcr_id,
            'file_type': file_type,
            'http_url': url,
            'ftp_url': ftp_url if ftp_url else '',
            'expected_filename': Path(path).name,
            'expected_path': str(path)
        })
    
    df_failed = pd.DataFrame(rows)
    df_failed.to_csv(failed_file, sep='\t', index=False)
    
    print(f"\n⚠ Created failed downloads list: {failed_file}")
    print(f"  Total failed files: {len(df_failed)}")
    print(f"  Unique failed runs: {df_failed['run_accession'].nunique()}")
    
    return failed_file

def main():
    parser = argparse.ArgumentParser(
        description='Download FASTQ files from ENA and create samples.tsv',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('ena_file', help='ENA file report (TSV format)')
    parser.add_argument('output_dir', help='Output directory for FASTQ files')
    parser.add_argument('--samples-tsv', default=None,
                       help='Output samples.tsv file path (default: <output_dir>/samples.tsv)')
    parser.add_argument('--dry-run', action='store_true',
                       help='Show what would be downloaded without actually downloading')
    parser.add_argument('--max-workers', type=int, default=8,
                       help='Maximum number of parallel downloads (default: 8)')
    parser.add_argument('--connections-per-file', type=int, default=2,
                       help='Number of connections per file for aria2c (default: 2, reduce if getting errors)')
    parser.add_argument('--no-aria2c', action='store_true',
                       help='Disable aria2c and use wget/curl with multiprocessing instead')
    
    args = parser.parse_args()
    
    # Check if ENA file exists
    if not os.path.exists(args.ena_file):
        print(f"ERROR: ENA file not found: {args.ena_file}")
        sys.exit(1)
    
    # Set default samples.tsv path
    if args.samples_tsv is None:
        args.samples_tsv = os.path.join(args.output_dir, 'samples.tsv')
    
    print("=" * 80)
    print("ENA FASTQ Downloader and samples.tsv Generator")
    print("=" * 80)
    print(f"ENA file: {args.ena_file}")
    print(f"Output directory: {args.output_dir}")
    print(f"Samples.tsv output: {args.samples_tsv}")
    print("=" * 80)
    
    # Parse ENA file
    print("\n[1/3] Parsing ENA file...")
    samples_data = parse_ena_file(args.ena_file)
    print(f"  Found {len(samples_data)} samples with valid FTP URLs")
    
    if len(samples_data) == 0:
        print("ERROR: No valid samples found in ENA file!")
        sys.exit(1)
    
    # Download files
    print("\n[2/3] Downloading FASTQ files...")
    downloaded_files, failed_runs = download_fastq_files(
        samples_data, 
        args.output_dir, 
        dry_run=args.dry_run,
        max_workers=args.max_workers,
        use_aria2c=not args.no_aria2c,
        connections_per_file=args.connections_per_file
    )
    
    if len(downloaded_files) == 0:
        print("ERROR: No files were downloaded!")
        sys.exit(1)
    
    # Write failed runs list if any
    if failed_runs:
        write_failed_runs(failed_runs, args.output_dir)
    
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

