#!/usr/bin/env python3
"""
Interactive Multi-Method Pathogen Detection Dashboard
Creates comprehensive interactive visualizations combining all detection methods
with filtering, sorting, confidence scores, and export capabilities.
"""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import plotly.offline as pyo
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import os
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

class PathogenDetectionDashboard:
    def __init__(self, results_dir="results"):
        self.results_dir = results_dir
        self.data = {}
        self.confidence_data = {}
        self.load_all_data()
    
    def load_all_data(self):
        """Load all detection method data"""
        print("Loading detection data...")
        
        # Load E-Score data
        self.load_escore_data()
        
        # Load HOPS data
        self.load_hops_data()
        
        # Load BWA data
        self.load_bwa_data()
        
        # Load comparison data
        self.load_comparison_data()
        
        # Load detection scores
        self.load_detection_scores()
        
        print(f"Loaded data for {len(self.data)} sample-pathogen pairs")
    
    def load_escore_data(self):
        """Load E-Score results"""
        escore_files = Path(self.results_dir).glob("*/evalue/pathogen/*_pathogen.csv")
        for file in escore_files:
            sample = file.parent.parent.name
            df = pd.read_csv(file)
            for _, row in df.iterrows():
                pathogen = row['taxonomy'].strip()
                key = f"{sample}_{pathogen}"
                self.data[key] = {
                    'sample': sample,
                    'pathogen': pathogen,
                    'escore': row.get('escore', 0),
                    'reads': row.get('reads', 0),
                    'coverage': row.get('coverage', 0)
                }
                # Calculate confidence based on E-Score
                self.confidence_data[key] = {
                    'escore_confidence': min(row.get('escore', 0) * 10, 1.0)
                }
    
    def load_hops_data(self):
        """Load HOPS results"""
        hops_file = Path(self.results_dir) / "hops/maltExtract/heatmap_overview_Wevid.tsv"
        if hops_file.exists():
            df = pd.read_csv(hops_file, sep='\t')
            df.columns = [col.strip().strip('"') for col in df.columns]
            
            for _, row in df.iterrows():
                pathogen = row['node'].strip('"')
                for col in df.columns[1:]:  # Skip 'node' column
                    if '_unaligned.rma6' in col:
                        sample = col.replace('_unaligned.rma6', '')
                        score = row[col]
                        key = f"{sample}_{pathogen}"
                        
                        if key in self.data:
                            self.data[key]['hops_score'] = score
                            # HOPS confidence based on score thresholds
                            if score >= 4:
                                self.confidence_data[key]['hops_confidence'] = 1.0
                            elif score >= 3:
                                self.confidence_data[key]['hops_confidence'] = 0.8
                            elif score >= 2:
                                self.confidence_data[key]['hops_confidence'] = 0.6
                            else:
                                self.confidence_data[key]['hops_confidence'] = 0.2
    
    def load_bwa_data(self):
        """Load BWA alignment data"""
        bwa_files = Path(self.results_dir).glob("*/pathogen_mapping/*.ani.txt")
        for file in bwa_files:
            parts = file.stem.split('_')
            if len(parts) >= 3:
                sample = parts[0]
                pathogen = '_'.join(parts[1:-1])  # Handle pathogens with underscores
                
                try:
                    with open(file, 'r') as f:
                        content = f.read().strip()
                        if "ANI ≈" in content:
                            ani = float(content.split("ANI ≈ ")[1].split("%")[0])
                            key = f"{sample}_{pathogen}"
                            
                            if key in self.data:
                                self.data[key]['ani'] = ani
                                # ANI confidence
                                if ani > 99:
                                    self.confidence_data[key]['ani_confidence'] = 1.0
                                elif ani > 97:
                                    self.confidence_data[key]['ani_confidence'] = 0.9
                                elif ani > 95:
                                    self.confidence_data[key]['ani_confidence'] = 0.7
                                else:
                                    self.confidence_data[key]['ani_confidence'] = 0.5
                except:
                    pass
    
    def load_comparison_data(self):
        """Load comparison results"""
        comparison_files = Path(self.results_dir).glob("comparison/*_comparison.tsv")
        for file in comparison_files:
            sample = file.stem.replace('_comparison', '')
            df = pd.read_csv(file, sep='\t')
            
            for _, row in df.iterrows():
                pathogen = row['Krakenuniq name']
                key = f"{sample}_{pathogen}"
                
                if key in self.data:
                    self.data[key]['in_krakenuniq'] = row.get('In_Krakenuniq', False)
                    self.data[key]['in_hops'] = row.get('In_HOPS', False)
    
    def load_detection_scores(self):
        """Load detection scores"""
        scores_file = Path(self.results_dir) / "pathogen_detection/detection_scores_matrix.csv"
        if scores_file.exists():
            df = pd.read_csv(scores_file, index_col=0)
            
            for sample in df.index:
                for pathogen in df.columns:
                    key = f"{sample}_{pathogen}"
                    if key in self.data:
                        self.data[key]['detection_score'] = df.loc[sample, pathogen]
    
    def create_interactive_dashboard(self):
        """Create the main interactive dashboard"""
        # Convert data to DataFrame
        df_data = []
        for key, values in self.data.items():
            row = values.copy()
            row['key'] = key
            
            # Add confidence scores
            if key in self.confidence_data:
                row.update(self.confidence_data[key])
            
            df_data.append(row)
        
        df = pd.DataFrame(df_data)
        
        if df.empty:
            return self.create_empty_dashboard()
        
        # Create subplots
        fig = make_subplots(
            rows=3, cols=2,
            subplot_titles=[
                'Detection Methods Comparison',
                'Confidence Scores Distribution',
                'E-Score vs HOPS Score',
                'ANI vs Detection Score',
                'Sample-Pathogen Heatmap',
                'Method Performance Summary'
            ],
            specs=[
                [{"type": "bar"}, {"type": "box"}],
                [{"type": "scatter"}, {"type": "scatter"}],
                [{"type": "heatmap"}, {"type": "bar"}]
            ],
            vertical_spacing=0.08,
            horizontal_spacing=0.1
        )
        
        # 1. Detection Methods Comparison
        methods_data = []
        for method in ['escore', 'hops_score', 'ani', 'detection_score']:
            if method in df.columns:
                methods_data.append(df[method].dropna())
        
        for i, (method, data) in enumerate(zip(['E-Score', 'HOPS', 'ANI', 'Detection Score'], methods_data)):
            fig.add_trace(
                go.Bar(
                    x=[method],
                    y=[data.mean()],
                    error_y=dict(type='data', array=[data.std()]),
                    name=method,
                    showlegend=False
                ),
                row=1, col=1
            )
        
        # 2. Confidence Scores Distribution
        confidence_cols = [col for col in df.columns if 'confidence' in col]
        for col in confidence_cols:
            fig.add_trace(
                go.Box(
                    y=df[col].dropna(),
                    name=col.replace('_confidence', ''),
                    showlegend=False
                ),
                row=1, col=2
            )
        
        # 3. E-Score vs HOPS Score
        if 'escore' in df.columns and 'hops_score' in df.columns:
            fig.add_trace(
                go.Scatter(
                    x=df['escore'],
                    y=df['hops_score'],
                    mode='markers',
                    text=df['pathogen'],
                    hovertemplate='<b>%{text}</b><br>E-Score: %{x}<br>HOPS: %{y}<extra></extra>',
                    name='Sample-Pathogen',
                    showlegend=False
                ),
                row=2, col=1
            )
        
        # 4. ANI vs Detection Score
        if 'ani' in df.columns and 'detection_score' in df.columns:
            fig.add_trace(
                go.Scatter(
                    x=df['ani'],
                    y=df['detection_score'],
                    mode='markers',
                    text=df['pathogen'],
                    hovertemplate='<b>%{text}</b><br>ANI: %{x}%<br>Score: %{y}<extra></extra>',
                    name='ANI-Score',
                    showlegend=False
                ),
                row=2, col=2
            )
        
        # 5. Sample-Pathogen Heatmap
        if 'detection_score' in df.columns:
            pivot_df = df.pivot_table(
                index='sample', 
                columns='pathogen', 
                values='detection_score', 
                fill_value=0
            )
            
            fig.add_trace(
                go.Heatmap(
                    z=pivot_df.values,
                    x=pivot_df.columns,
                    y=pivot_df.index,
                    colorscale='RdYlBu_r',
                    showscale=True,
                    name='Detection Scores'
                ),
                row=3, col=1
            )
        
        # 6. Method Performance Summary
        method_performance = []
        for method in ['escore', 'hops_score', 'ani', 'detection_score']:
            if method in df.columns:
                method_performance.append({
                    'method': method.replace('_', ' ').title(),
                    'mean': df[method].mean(),
                    'std': df[method].std(),
                    'count': df[method].count()
                })
        
        perf_df = pd.DataFrame(method_performance)
        fig.add_trace(
            go.Bar(
                x=perf_df['method'],
                y=perf_df['mean'],
                error_y=dict(type='data', array=perf_df['std']),
                text=perf_df['count'],
                textposition='auto',
                hovertemplate='<b>%{x}</b><br>Mean: %{y:.2f}<br>Samples: %{text}<extra></extra>',
                showlegend=False
            ),
            row=3, col=2
        )
        
        # Update layout
        fig.update_layout(
            title={
                'text': 'PIGSTI Pathogen Detection Interactive Dashboard',
                'x': 0.5,
                'xanchor': 'center',
                'font': {'size': 20}
            },
            height=1200,
            showlegend=False,
            template='plotly_white'
        )
        
        # Update axes labels
        fig.update_xaxes(title_text="Detection Methods", row=1, col=1)
        fig.update_yaxes(title_text="Score", row=1, col=1)
        
        fig.update_yaxes(title_text="Confidence Score", row=1, col=2)
        
        fig.update_xaxes(title_text="E-Score", row=2, col=1)
        fig.update_yaxes(title_text="HOPS Score", row=2, col=1)
        
        fig.update_xaxes(title_text="ANI (%)", row=2, col=2)
        fig.update_yaxes(title_text="Detection Score", row=2, col=2)
        
        fig.update_xaxes(title_text="Pathogen", row=3, col=1)
        fig.update_yaxes(title_text="Sample", row=3, col=1)
        
        fig.update_xaxes(title_text="Method", row=3, col=2)
        fig.update_yaxes(title_text="Mean Score", row=3, col=2)
        
        return fig
    
    def create_empty_dashboard(self):
        """Create empty dashboard when no data is available"""
        fig = go.Figure()
        fig.add_annotation(
            text="No pathogen detection data available",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=20)
        )
        fig.update_layout(
            title="PIGSTI Pathogen Detection Dashboard",
            xaxis=dict(showgrid=False, showticklabels=False),
            yaxis=dict(showgrid=False, showticklabels=False)
        )
        return fig
    
    def create_filterable_table(self):
        """Create filterable data table"""
        df_data = []
        for key, values in self.data.items():
            row = values.copy()
            row['key'] = key
            
            # Add confidence scores
            if key in self.confidence_data:
                row.update(self.confidence_data[key])
            
            df_data.append(row)
        
        df = pd.DataFrame(df_data)
        
        if df.empty:
            return go.Figure()
        
        # Create interactive table
        fig = go.Figure(data=[go.Table(
            header=dict(
                values=list(df.columns),
                fill_color='paleturquoise',
                align='left',
                font=dict(size=12, color='black')
            ),
            cells=dict(
                values=[df[col] for col in df.columns],
                fill_color='lavender',
                align='left',
                font=dict(size=11)
            )
        )])
        
        fig.update_layout(
            title="Pathogen Detection Results - Filterable Table",
            height=600
        )
        
        return fig
    
    def generate_dashboard_html(self, output_file="results/final/pathogen_detection_dashboard.html"):
        """Generate the complete interactive dashboard HTML"""
        print("Generating interactive dashboard...")
        
        # Create main dashboard
        main_fig = self.create_interactive_dashboard()
        
        # Create filterable table
        table_fig = self.create_filterable_table()
        
        # Generate HTML content
        html_content = f"""
        <!DOCTYPE html>
        <html>
        <head>
            <title>PIGSTI Pathogen Detection Dashboard</title>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                .header {{ background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                          color: white; padding: 20px; border-radius: 10px; margin-bottom: 20px; }}
                .section {{ margin: 30px 0; }}
                .stats {{ display: flex; justify-content: space-around; margin: 20px 0; }}
                .stat-box {{ background: #f8f9fa; padding: 15px; border-radius: 8px; text-align: center; }}
                .controls {{ background: #e9ecef; padding: 15px; border-radius: 8px; margin: 20px 0; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>🧬 PIGSTI Pathogen Detection Dashboard</h1>
                <p>Interactive Multi-Method Pathogen Detection Analysis</p>
                <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            </div>
            
            <div class="stats">
                <div class="stat-box">
                    <h3>{len(self.data)}</h3>
                    <p>Sample-Pathogen Pairs</p>
                </div>
                <div class="stat-box">
                    <h3>{len(set([d['sample'] for d in self.data.values()]))}</h3>
                    <p>Samples Analyzed</p>
                </div>
                <div class="stat-box">
                    <h3>{len(set([d['pathogen'] for d in self.data.values()]))}</h3>
                    <p>Pathogens Detected</p>
                </div>
            </div>
            
            <div class="section">
                <h2>📊 Interactive Analysis Dashboard</h2>
                {main_fig.to_html(full_html=False, include_plotlyjs='cdn')}
            </div>
            
            <div class="section">
                <h2>📋 Detailed Results Table</h2>
                {table_fig.to_html(full_html=False, include_plotlyjs='cdn')}
            </div>
            
            <div class="controls">
                <h3>🎛️ Export Options</h3>
                <p>Use the toolbar above each plot to:</p>
                <ul>
                    <li><strong>📥 Download as PNG/SVG/PDF</strong> - High-resolution images for presentations</li>
                    <li><strong>📊 Download as HTML</strong> - Interactive plots for web sharing</li>
                    <li><strong>🔍 Zoom and Pan</strong> - Explore data in detail</li>
                    <li><strong>📱 Responsive View</strong> - Optimized for different screen sizes</li>
                </ul>
            </div>
        </body>
        </html>
        """
        
        # Save HTML file
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        with open(output_file, 'w') as f:
            f.write(html_content)
        
        print(f"Interactive dashboard saved to: {output_file}")
        return output_file

def main():
    """Main function to generate the dashboard"""
    dashboard = PathogenDetectionDashboard()
    dashboard.generate_dashboard_html()

if __name__ == "__main__":
    main()
