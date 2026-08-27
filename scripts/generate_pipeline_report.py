#!/usr/bin/env python3
"""
Pipeline Execution Report Generator
Generates comprehensive reports including:
- Execution timing and status heatmap
- Pipeline workflow diagram
- Step completion analysis
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime
import re
import glob
from pathlib import Path
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.offline as pyo
import networkx as nx
import warnings
warnings.filterwarnings('ignore')

def parse_log_timing(log_file):
    """Parse log file to extract start and end times"""
    try:
        with open(log_file, 'r') as f:
            content = f.read()
        
        # Extract start and end times
        start_match = re.search(r'Starting.*at (.*)', content)
        end_match = re.search(r'completed.*at (.*)', content)
        
        if start_match and end_match:
            start_time = datetime.strptime(start_match.group(1).strip(), '%a %b %d %H:%M:%S %Z %Y')
            end_time = datetime.strptime(end_match.group(1).strip(), '%a %b %d %H:%M:%S %Z %Y')
            duration = (end_time - start_time).total_seconds() / 60  # minutes
            return start_time, end_time, duration
    except:
        pass
    return None, None, None

def check_file_completion(file_path):
    """Check if a file exists and get its modification time"""
    if os.path.exists(file_path):
        mtime = os.path.getmtime(file_path)
        return True, datetime.fromtimestamp(mtime)
    return False, None

def generate_timing_data():
    """Generate timing data from log files"""
    timing_data = []
    
    # Get all log files (only if they exist)
    log_patterns = [
        'logs/adapter_removal/*.log',
        'logs/krakenuniq/*.log',
        'logs/evalue/*.log',
        'logs/host_mapping/*.log',
        'logs/mtdna_mapping/*.log',
        'logs/pathogen_mapping/*.log'
    ]
    
    for pattern in log_patterns:
        log_files = glob.glob(pattern)
        if not log_files:
            print(f"No log files found for pattern: {pattern}")
            continue
            
        for log_file in log_files:
            if os.path.exists(log_file):
                start_time, end_time, duration = parse_log_timing(log_file)
                if start_time and end_time:
                    # Extract sample and step from filename
                    filename = os.path.basename(log_file)
                    if '_pe.log' in filename or '_se.log' in filename:
                        sample = filename.replace('_pe.log', '').replace('_se.log', '')
                        step = 'adapter_removal'
                    else:
                        sample = filename.replace('.log', '')
                        step = os.path.dirname(log_file).split('/')[-1]
                    
                    timing_data.append({
                        'sample': sample,
                        'step': step,
                        'start_time': start_time,
                        'end_time': end_time,
                        'duration_minutes': duration,
                        'log_file': log_file
                    })
    
    return pd.DataFrame(timing_data)

def generate_completion_matrix():
    """Generate completion matrix for all samples and steps"""
    # Define expected steps and their output files
    steps = {
        'adapter_removal': 'results/libraries/{sample}/adapter_removal/{sample}.collapsed.gz',
        'fastq_screen': 'results/libraries/{sample}/fastq_screen/{sample}_best_species.txt',
        'krakenuniq': 'results/metagenomics/krakenuniq/{sample}/{sample}_kraken-report.txt',
        'escore': 'results/pathogen/{sample}/evalue/pathogen/{sample}_pathogen.csv',
        'host_alignment': 'results/libraries/{sample}/host_mapping/{sample}.dedup.bam',
        'mtdna_alignment': 'results/libraries/{sample}/mtdna_mapping/{sample}.dedup.bam',
        'pathogen_alignment': 'results/pathogen/{sample}/pathogen_mapping/',
        'damageprofiler_host': 'results/libraries/{sample}/damageprofiler_host/',
        'damageprofiler_mtdna': 'results/libraries/{sample}/damageprofiler_mtdna/',
        'qualimap_host': 'results/libraries/{sample}/qualimap/genome_results.txt',
        'qualimap_mtdna': 'results/libraries/{sample}/qualimap_mtdna/genome_results.txt'
    }
    
    # Get all samples from results directory (sorted for stable reports)
    samples = []
    results_dir = 'results'
    if not os.path.isdir(results_dir):
        return pd.DataFrame()
    for item in sorted(os.listdir(results_dir)):
        if os.path.isdir(os.path.join(results_dir, item)) and not item.startswith('.'):
            samples.append(item)
    
    # Create completion matrix
    completion_matrix = []
    for sample in samples:
        row = {'sample': sample}
        for step, pattern in steps.items():
            if '{sample}' in pattern:
                file_path = pattern.format(sample=sample)
            else:
                file_path = pattern
            
            completed, mtime = check_file_completion(file_path)
            row[step] = 1 if completed else 0
            row[f'{step}_time'] = mtime
        
        completion_matrix.append(row)
    
    return pd.DataFrame(completion_matrix)

def create_execution_heatmap(timing_df, completion_df):
    """Create execution heatmap showing timing and completion status"""
    if timing_df.empty:
        print("No timing data available")
        return
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=['Step Completion Status', 'Execution Duration (minutes)'],
        vertical_spacing=0.1
    )
    
    # Completion heatmap
    completion_matrix = completion_df.set_index('sample').drop(columns=[col for col in completion_df.columns if '_time' in col])
    
    fig.add_trace(
        go.Heatmap(
            z=completion_matrix.values,
            x=completion_matrix.columns,
            y=completion_matrix.index,
            colorscale='RdYlGn',
            showscale=True,
            name='Completion Status'
        ),
        row=1, col=1
    )
    
    # Duration heatmap
    if not timing_df.empty:
        pivot_duration = timing_df.pivot_table(
            index='sample', 
            columns='step', 
            values='duration_minutes', 
            fill_value=0
        )
        
        fig.add_trace(
            go.Heatmap(
                z=pivot_duration.values,
                x=pivot_duration.columns,
                y=pivot_duration.index,
                colorscale='Blues',
                showscale=True,
                name='Duration (min)'
            ),
            row=2, col=1
        )
    
    fig.update_layout(
        title='PIGSTI Pipeline Execution Report',
        height=800,
        showlegend=False
    )
    
    return fig

def create_workflow_diagram():
    """
    Build a static workflow figure (deterministic layout).
    Delegates to scripts/render_workflow_diagram.py — no random graph layout.
    """
    import sys

    scripts_dir = Path(__file__).resolve().parent
    if str(scripts_dir) not in sys.path:
        sys.path.insert(0, str(scripts_dir))
    from render_workflow_diagram import render

    png_path = Path("results/final/pipeline_workflow_diagram.png")
    svg_path = Path("results/final/pipeline_workflow_diagram.svg")
    render(png_path, svg_path)

    # Embed PNG in HTML report via a simple Plotly image trace
    import base64

    b64 = base64.b64encode(png_path.read_bytes()).decode("ascii")
    fig = go.Figure()
    fig.add_layout_image(
        dict(
            source=f"data:image/png;base64,{b64}",
            xref="paper",
            yref="paper",
            x=0,
            y=1,
            sizex=1,
            sizey=1,
            sizing="contain",
            layer="below",
        )
    )
    fig.update_layout(
        title="PIGSTI Pipeline Workflow",
        xaxis=dict(visible=False, range=[0, 1]),
        yaxis=dict(visible=False, range=[0, 1]),
        margin=dict(l=10, r=10, t=40, b=10),
        height=700,
    )
    return fig

def generate_html_report(timing_df, completion_df, heatmap_fig, workflow_fig):
    """Generate comprehensive HTML report"""
    
    # Calculate summary statistics
    total_samples = len(completion_df)
    
    # Only sum numeric columns (exclude datetime columns)
    numeric_cols = completion_df.select_dtypes(include=[np.number]).columns
    if 'sample' in numeric_cols:
        numeric_cols = numeric_cols.drop('sample')
    
    if len(numeric_cols) > 0:
        completed_steps = completion_df[numeric_cols].sum().sum()
        total_possible_steps = total_samples * len(numeric_cols)
    else:
        completed_steps = 0
        total_possible_steps = 0
    
    completion_rate = (completed_steps / total_possible_steps) * 100 if total_possible_steps > 0 else 0
    
    # Generate HTML content
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>PIGSTI Pipeline Execution Report</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 40px; }}
            .header {{ background-color: #f0f0f0; padding: 20px; border-radius: 5px; }}
            .summary {{ background-color: #e8f4fd; padding: 15px; margin: 20px 0; border-radius: 5px; }}
            .section {{ margin: 30px 0; }}
            table {{ border-collapse: collapse; width: 100%; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
            th {{ background-color: #f2f2f2; }}
        </style>
    </head>
    <body>
        <div class="header">
            <h1>PIGSTI Pipeline Execution Report</h1>
            <p>Generated on: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <div class="summary">
            <h2>Pipeline Summary</h2>
            <ul>
                <li><strong>Total Samples:</strong> {total_samples}</li>
                <li><strong>Completion Rate:</strong> {completion_rate:.1f}%</li>
                <li><strong>Completed Steps:</strong> {completed_steps} / {total_possible_steps}</li>
            </ul>
        </div>
        
        <div class="section">
            <h2>Execution Heatmap</h2>
            {heatmap_fig.to_html(full_html=False, include_plotlyjs='cdn') if heatmap_fig is not None else '<p>No timing data available for heatmap</p>'}
        </div>
        
        <div class="section">
            <h2>Pipeline Workflow Diagram</h2>
            {workflow_fig.to_html(full_html=False, include_plotlyjs='cdn') if workflow_fig is not None else '<p>Workflow diagram not available</p>'}
        </div>
        
        <div class="section">
            <h2>Sample Completion Status</h2>
            {completion_df.to_html(index=False, classes='table')}
        </div>
        
        <div class="section">
            <h2>Timing Data</h2>
            {timing_df.to_html(index=False, classes='table') if not timing_df.empty else '<p>No timing data available</p>'}
        </div>
    </body>
    </html>
    """
    
    return html_content

def main():
    """Main function to generate pipeline report"""
    print("Generating pipeline execution report...")
    
    # Create output directories
    os.makedirs('results', exist_ok=True)
    os.makedirs('logs', exist_ok=True)
    
    # Generate data
    print("Collecting timing data...")
    timing_df = generate_timing_data()
    
    print("Generating completion matrix...")
    completion_df = generate_completion_matrix()
    
    # Save timing data (even if empty)
    timing_df.to_csv('results/final/pipeline_timing_data.csv', index=False)
    if timing_df.empty:
        print("No timing data available - log files not found yet")
    else:
        print(f"Timing data saved to results/final/pipeline_timing_data.csv ({len(timing_df)} entries)")
    
    # Create visualizations
    print("Creating execution heatmap...")
    heatmap_fig = create_execution_heatmap(timing_df, completion_df)
    
    print("Creating workflow diagram...")
    workflow_fig = create_workflow_diagram()
    workflow_fig.write_html('results/final/pipeline_workflow_diagram.html')
    print("Workflow diagram saved to results/final/pipeline_workflow_diagram.png")
    print("Workflow diagram saved to results/final/pipeline_workflow_diagram.svg")
    print("Workflow diagram saved as HTML: results/final/pipeline_workflow_diagram.html")
    # Keep docs/ in sync for README (publication figure)
    docs_png = Path("docs/images/pipeline_workflow.png")
    docs_svg = Path("docs/images/pipeline_workflow.svg")
    if Path("results/final/pipeline_workflow_diagram.png").exists():
        docs_png.parent.mkdir(parents=True, exist_ok=True)
        docs_png.write_bytes(Path("results/final/pipeline_workflow_diagram.png").read_bytes())
        if Path("results/final/pipeline_workflow_diagram.svg").exists():
            docs_svg.write_text(
                Path("results/final/pipeline_workflow_diagram.svg").read_text(encoding="utf-8"),
                encoding="utf-8",
            )
    
    # Generate HTML report
    print("Generating HTML report...")
    html_content = generate_html_report(timing_df, completion_df, heatmap_fig, workflow_fig)
    
    with open('results/final/pipeline_execution_report.html', 'w') as f:
        f.write(html_content)
    
    print("Pipeline execution report generated successfully!")
    print("Files created:")
    print("- results/final/pipeline_execution_report.html")
    print("- results/final/pipeline_timing_data.csv")
    print("- results/final/pipeline_workflow_diagram.png")

if __name__ == "__main__":
    main()
