#!/usr/bin/env python3
"""
Publication-Ready Output Suite
Creates high-quality, publication-ready figures with consistent styling,
multi-panel layouts, and automatic caption generation.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
from pathlib import Path
import os
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Set publication-quality style
plt.style.use('seaborn-v0_8-whitegrid')
sns.set_palette("husl")

class PublicationFigureGenerator:
    def __init__(self, results_dir="results", output_dir="results/publication_figures"):
        self.results_dir = results_dir
        self.output_dir = output_dir
        self.data = {}
        self.load_all_data()
        
        # Publication settings
        self.dpi = 300
        self.font_size = 12
        self.title_size = 16
        self.label_size = 14
        
        # Color schemes
        self.colors = {
            'primary': '#2E86AB',
            'secondary': '#A23B72', 
            'accent': '#F18F01',
            'success': '#C73E1D',
            'neutral': '#6C757D'
        }
        
        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)
    
    def load_all_data(self):
        """Load all detection method data"""
        print("Loading data for publication figures...")
        
        # Load detection scores
        scores_file = Path(self.results_dir) / "pathogen_detection/detection_scores_matrix.csv"
        if scores_file.exists():
            self.data['scores'] = pd.read_csv(scores_file, index_col=0)
        
        # Load detailed scores
        detailed_file = Path(self.results_dir) / "pathogen_detection/detailed_scores.csv"
        if detailed_file.exists():
            self.data['detailed'] = pd.read_csv(detailed_file)
        
        # Load comparison data
        self.load_comparison_data()
        
        # Load abundance matrices
        self.load_abundance_data()
    
    def load_comparison_data(self):
        """Load comparison data from all samples"""
        comparison_files = Path(self.results_dir).glob("comparison/*_comparison.tsv")
        comparison_data = []
        
        for file in comparison_files:
            sample = file.stem.replace('_comparison', '')
            df = pd.read_csv(file, sep='\t')
            df['sample'] = sample
            comparison_data.append(df)
        
        if comparison_data:
            self.data['comparison'] = pd.concat(comparison_data, ignore_index=True)
    
    def load_abundance_data(self):
        """Load abundance matrix data"""
        abundance_file = Path(self.results_dir) / "KRAKENUNIQ_ABUNDANCE_MATRIX/krakenuniq_abundance_matrix_absolute.csv"
        if abundance_file.exists():
            self.data['abundance'] = pd.read_csv(abundance_file, index_col=0)
    
    def create_multi_method_comparison(self):
        """Create multi-panel figure comparing all detection methods"""
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Multi-Method Pathogen Detection Analysis', fontsize=self.title_size, fontweight='bold')
        
        # Panel A: Detection Score Distribution
        if 'scores' in self.data:
            scores_flat = self.data['scores'].values.flatten()
            scores_flat = scores_flat[scores_flat > 0]  # Remove zeros
            
            axes[0, 0].hist(scores_flat, bins=20, alpha=0.7, color=self.colors['primary'], edgecolor='black')
            axes[0, 0].set_xlabel('Detection Score', fontsize=self.label_size)
            axes[0, 0].set_ylabel('Frequency', fontsize=self.label_size)
            axes[0, 0].set_title('A) Detection Score Distribution', fontsize=self.label_size, fontweight='bold')
            axes[0, 0].grid(True, alpha=0.3)
        
        # Panel B: Method Performance Comparison
        if 'detailed' in self.data:
            method_cols = [col for col in self.data['detailed'].columns if col not in ['sample', 'pathogen']]
            method_performance = self.data['detailed'][method_cols].mean()
            
            bars = axes[0, 1].bar(range(len(method_performance)), method_performance.values, 
                                color=[self.colors['primary'], self.colors['secondary'], 
                                      self.colors['accent'], self.colors['success']][:len(method_performance)])
            axes[0, 1].set_xlabel('Detection Criteria', fontsize=self.label_size)
            axes[0, 1].set_ylabel('Mean Score', fontsize=self.label_size)
            axes[0, 1].set_title('B) Method Performance Comparison', fontsize=self.label_size, fontweight='bold')
            axes[0, 1].set_xticks(range(len(method_performance)))
            axes[0, 1].set_xticklabels([col.replace('_', ' ').title() for col in method_performance.index], 
                                     rotation=45, ha='right')
            axes[0, 1].grid(True, alpha=0.3)
            
            # Add value labels on bars
            for bar, value in zip(bars, method_performance.values):
                axes[0, 1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                              f'{value:.2f}', ha='center', va='bottom', fontsize=10)
        
        # Panel C: Sample-Pathogen Heatmap
        if 'scores' in self.data:
            im = axes[1, 0].imshow(self.data['scores'].values, cmap='RdYlBu_r', aspect='auto')
            axes[1, 0].set_xlabel('Pathogen', fontsize=self.label_size)
            axes[1, 0].set_ylabel('Sample', fontsize=self.label_size)
            axes[1, 0].set_title('C) Detection Score Heatmap', fontsize=self.label_size, fontweight='bold')
            
            # Set tick labels
            axes[1, 0].set_xticks(range(len(self.data['scores'].columns)))
            axes[1, 0].set_xticklabels(self.data['scores'].columns, rotation=45, ha='right')
            axes[1, 0].set_yticks(range(len(self.data['scores'].index)))
            axes[1, 0].set_yticklabels(self.data['scores'].index)
            
            # Add colorbar
            cbar = plt.colorbar(im, ax=axes[1, 0])
            cbar.set_label('Detection Score', fontsize=self.label_size)
        
        # Panel D: Confidence Analysis
        if 'detailed' in self.data:
            # Calculate overall confidence for each sample-pathogen pair
            confidence_data = []
            for _, row in self.data['detailed'].iterrows():
                sample = row['sample']
                pathogen = row['pathogen']
                # Calculate confidence as mean of all criteria
                criteria_scores = [row[col] for col in method_cols if col in row]
                confidence = np.mean(criteria_scores) if criteria_scores else 0
                confidence_data.append({'sample': sample, 'pathogen': pathogen, 'confidence': confidence})
            
            conf_df = pd.DataFrame(confidence_data)
            if not conf_df.empty:
                conf_pivot = conf_df.pivot(index='sample', columns='pathogen', values='confidence')
                
                im = axes[1, 1].imshow(conf_pivot.values, cmap='viridis', aspect='auto')
                axes[1, 1].set_xlabel('Pathogen', fontsize=self.label_size)
                axes[1, 1].set_ylabel('Sample', fontsize=self.label_size)
                axes[1, 1].set_title('D) Detection Confidence Analysis', fontsize=self.label_size, fontweight='bold')
                
                # Set tick labels
                axes[1, 1].set_xticks(range(len(conf_pivot.columns)))
                axes[1, 1].set_xticklabels(conf_pivot.columns, rotation=45, ha='right')
                axes[1, 1].set_yticks(range(len(conf_pivot.index)))
                axes[1, 1].set_yticklabels(conf_pivot.index)
                
                # Add colorbar
                cbar = plt.colorbar(im, ax=axes[1, 1])
                cbar.set_label('Confidence Score', fontsize=self.label_size)
        
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(self.output_dir, "multi_method_comparison.png")
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.savefig(output_file.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
        plt.close()  # Close the figure to free memory
        
        # Generate caption
        caption = self.generate_caption("multi_method_comparison", {
            'total_samples': len(self.data['scores'].index) if 'scores' in self.data else 0,
            'total_pathogens': len(self.data['scores'].columns) if 'scores' in self.data else 0,
            'detection_pairs': np.sum(self.data['scores'].values > 0) if 'scores' in self.data else 0
        })
        
        self.save_caption("multi_method_comparison", caption)
        
        print(f"Multi-method comparison figure saved to: {output_file}")
        return output_file
    
    def create_method_performance_analysis(self):
        """Create detailed method performance analysis"""
        if 'detailed' not in self.data:
            print("No detailed scores data available")
            return None
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Detection Method Performance Analysis', fontsize=self.title_size, fontweight='bold')
        
        method_cols = [col for col in self.data['detailed'].columns if col not in ['sample', 'pathogen']]
        
        # Panel A: Individual Method Scores
        for i, method in enumerate(method_cols[:4]):  # Limit to 4 methods
            row, col = i // 2, i % 2
            scores = self.data['detailed'][method].dropna()
            
            axes[row, col].hist(scores, bins=15, alpha=0.7, color=self.colors['primary'], edgecolor='black')
            axes[row, col].set_xlabel('Score', fontsize=self.label_size)
            axes[row, col].set_ylabel('Frequency', fontsize=self.label_size)
            axes[row, col].set_title(f'{chr(65+i)}) {method.replace("_", " ").title()}', 
                                   fontsize=self.label_size, fontweight='bold')
            axes[row, col].grid(True, alpha=0.3)
            
            # Add statistics
            mean_score = scores.mean()
            axes[row, col].axvline(mean_score, color='red', linestyle='--', linewidth=2, 
                                 label=f'Mean: {mean_score:.2f}')
            axes[row, col].legend()
        
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(self.output_dir, "method_performance_analysis.png")
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.savefig(output_file.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
        plt.close()  # Close the figure to free memory
        
        # Generate caption
        caption = self.generate_caption("method_performance_analysis", {
            'methods_analyzed': len(method_cols),
            'total_assessments': len(self.data['detailed'])
        })
        
        self.save_caption("method_performance_analysis", caption)
        
        print(f"Method performance analysis saved to: {output_file}")
        return output_file
    
    def create_sample_comparison(self):
        """Create sample comparison analysis"""
        if 'scores' not in self.data:
            print("No scores data available")
            return None
        
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('Sample-Level Pathogen Detection Analysis', fontsize=self.title_size, fontweight='bold')
        
        # Panel A: Sample Detection Counts
        sample_counts = (self.data['scores'] > 0).sum(axis=1)
        axes[0, 0].bar(range(len(sample_counts)), sample_counts.values, 
                      color=self.colors['primary'], alpha=0.7)
        axes[0, 0].set_xlabel('Sample', fontsize=self.label_size)
        axes[0, 0].set_ylabel('Pathogens Detected', fontsize=self.label_size)
        axes[0, 0].set_title('A) Pathogens Detected per Sample', fontsize=self.label_size, fontweight='bold')
        axes[0, 0].set_xticks(range(len(sample_counts)))
        axes[0, 0].set_xticklabels(sample_counts.index, rotation=45, ha='right')
        axes[0, 0].grid(True, alpha=0.3)
        
        # Panel B: Sample Detection Scores Distribution
        sample_means = self.data['scores'].mean(axis=1)
        axes[0, 1].hist(sample_means, bins=10, alpha=0.7, color=self.colors['secondary'], edgecolor='black')
        axes[0, 1].set_xlabel('Mean Detection Score', fontsize=self.label_size)
        axes[0, 1].set_ylabel('Frequency', fontsize=self.label_size)
        axes[0, 1].set_title('B) Sample Detection Score Distribution', fontsize=self.label_size, fontweight='bold')
        axes[0, 1].grid(True, alpha=0.3)
        
        # Panel C: Pathogen Prevalence
        pathogen_counts = (self.data['scores'] > 0).sum(axis=0)
        axes[1, 0].bar(range(len(pathogen_counts)), pathogen_counts.values, 
                      color=self.colors['accent'], alpha=0.7)
        axes[1, 0].set_xlabel('Pathogen', fontsize=self.label_size)
        axes[1, 0].set_ylabel('Samples with Detection', fontsize=self.label_size)
        axes[1, 0].set_title('C) Pathogen Prevalence Across Samples', fontsize=self.label_size, fontweight='bold')
        axes[1, 0].set_xticks(range(len(pathogen_counts)))
        axes[1, 0].set_xticklabels(pathogen_counts.index, rotation=45, ha='right')
        axes[1, 0].grid(True, alpha=0.3)
        
        # Panel D: Detection Score vs Prevalence
        axes[1, 1].scatter(pathogen_counts.values, self.data['scores'].mean(axis=0).values,
                          alpha=0.7, s=100, color=self.colors['success'])
        axes[1, 1].set_xlabel('Prevalence (Number of Samples)', fontsize=self.label_size)
        axes[1, 1].set_ylabel('Mean Detection Score', fontsize=self.label_size)
        axes[1, 1].set_title('D) Detection Score vs Prevalence', fontsize=self.label_size, fontweight='bold')
        axes[1, 1].grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(self.output_dir, "sample_comparison.png")
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.savefig(output_file.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
        plt.close()  # Close the figure to free memory
        
        # Generate caption
        caption = self.generate_caption("sample_comparison", {
            'total_samples': len(self.data['scores'].index),
            'total_pathogens': len(self.data['scores'].columns),
            'mean_pathogens_per_sample': sample_counts.mean()
        })
        
        self.save_caption("sample_comparison", caption)
        
        print(f"Sample comparison analysis saved to: {output_file}")
        return output_file
    
    def create_summary_statistics(self):
        """Create summary statistics figure"""
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle('PIGSTI Pipeline Summary Statistics', fontsize=self.title_size, fontweight='bold')
        
        # Panel A: Overall Detection Summary
        if 'scores' in self.data:
            total_detections = np.sum(self.data['scores'].values > 0)
            total_possible = self.data['scores'].size
            detection_rate = total_detections / total_possible * 100
            
            labels = ['Detected', 'Not Detected']
            sizes = [total_detections, total_possible - total_detections]
            colors = [self.colors['success'], self.colors['neutral']]
            
            axes[0, 0].pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
            axes[0, 0].set_title(f'A) Overall Detection Rate\n({detection_rate:.1f}%)', 
                               fontsize=self.label_size, fontweight='bold')
        
        # Panel B: Score Distribution
        if 'scores' in self.data:
            scores_flat = self.data['scores'].values.flatten()
            scores_flat = scores_flat[scores_flat > 0]
            
            axes[0, 1].hist(scores_flat, bins=20, alpha=0.7, color=self.colors['primary'], edgecolor='black')
            axes[0, 1].set_xlabel('Detection Score', fontsize=self.label_size)
            axes[0, 1].set_ylabel('Frequency', fontsize=self.label_size)
            axes[0, 1].set_title('B) Detection Score Distribution', fontsize=self.label_size, fontweight='bold')
            axes[0, 1].grid(True, alpha=0.3)
        
        # Panel C: Method Comparison (if detailed data available)
        if 'detailed' in self.data:
            method_cols = [col for col in self.data['detailed'].columns if col not in ['sample', 'pathogen']]
            method_means = self.data['detailed'][method_cols].mean()
            
            bars = axes[1, 0].bar(range(len(method_means)), method_means.values,
                                color=[self.colors['primary'], self.colors['secondary'], 
                                      self.colors['accent'], self.colors['success']][:len(method_means)])
            axes[1, 0].set_xlabel('Detection Method', fontsize=self.label_size)
            axes[1, 0].set_ylabel('Mean Score', fontsize=self.label_size)
            axes[1, 0].set_title('C) Method Performance Comparison', fontsize=self.label_size, fontweight='bold')
            axes[1, 0].set_xticks(range(len(method_means)))
            axes[1, 0].set_xticklabels([col.replace('_', ' ').title() for col in method_means.index], 
                                     rotation=45, ha='right')
            axes[1, 0].grid(True, alpha=0.3)
        
        # Panel D: Quality Metrics
        if 'scores' in self.data:
            high_confidence = np.sum(self.data['scores'].values >= 7)  # High confidence threshold
            medium_confidence = np.sum((self.data['scores'].values >= 4) & (self.data['scores'].values < 7))
            low_confidence = np.sum((self.data['scores'].values > 0) & (self.data['scores'].values < 4))
            
            categories = ['High\n(≥7)', 'Medium\n(4-6)', 'Low\n(1-3)']
            values = [high_confidence, medium_confidence, low_confidence]
            colors = [self.colors['success'], self.colors['accent'], self.colors['secondary']]
            
            bars = axes[1, 1].bar(categories, values, color=colors, alpha=0.7)
            axes[1, 1].set_ylabel('Number of Detections', fontsize=self.label_size)
            axes[1, 1].set_title('D) Detection Confidence Levels', fontsize=self.label_size, fontweight='bold')
            axes[1, 1].grid(True, alpha=0.3)
            
            # Add value labels
            for bar, value in zip(bars, values):
                axes[1, 1].text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                              str(value), ha='center', va='bottom', fontsize=12, fontweight='bold')
        
        plt.tight_layout()
        
        # Save figure
        output_file = os.path.join(self.output_dir, "summary_statistics.png")
        plt.savefig(output_file, dpi=self.dpi, bbox_inches='tight', facecolor='white')
        plt.savefig(output_file.replace('.png', '.pdf'), bbox_inches='tight', facecolor='white')
        plt.close()  # Close the figure to free memory
        
        # Generate caption
        caption = self.generate_caption("summary_statistics", {
            'total_samples': len(self.data['scores'].index) if 'scores' in self.data else 0,
            'total_pathogens': len(self.data['scores'].columns) if 'scores' in self.data else 0,
            'detection_rate': detection_rate if 'scores' in self.data else 0
        })
        
        self.save_caption("summary_statistics", caption)
        
        print(f"Summary statistics saved to: {output_file}")
        return output_file
    
    def generate_caption(self, figure_type, stats):
        """Generate automatic captions for figures"""
        templates = {
            "multi_method_comparison": (
                "\n            Figure 1. Multi-method pathogen detection analysis. (A) Distribution of detection scores across all sample-pathogen pairs (n={detection_pairs}).\n"
                "            (B) Performance comparison of different detection criteria. (C) Heatmap showing detection scores for each sample-pathogen combination\n"
                "            (n={total_samples} samples, {total_pathogens} pathogens). (D) Confidence analysis based on multiple detection criteria.\n"
                "            Higher scores indicate stronger evidence for pathogen presence.\n            "
            ),
            "method_performance_analysis": (
                "\n            Figure 2. Detection method performance analysis. Individual performance metrics for each detection criterion across\n"
                "            {total_assessments} assessments. Methods analyzed: {methods_analyzed} different criteria.\n"
                "            Higher scores indicate better performance for each criterion.\n            "
            ),
            "sample_comparison": (
                "\n            Figure 3. Sample-level pathogen detection analysis. (A) Number of pathogens detected per sample (n={total_samples} samples).\n"
                "            (B) Distribution of mean detection scores across samples. (C) Prevalence of each pathogen across samples (n={total_pathogens} pathogens).\n"
                "            (D) Relationship between pathogen prevalence and mean detection score. Mean pathogens per sample: {mean_pathogens_per_sample:.1f}.\n            "
            ),
            "summary_statistics": (
                "\n            Figure 4. PIGSTI pipeline summary statistics. (A) Overall detection rate across all sample-pathogen combinations\n"
                "            (n={total_samples} samples, {total_pathogens} pathogens). (B) Distribution of detection scores for positive detections.\n"
                "            (C) Performance comparison of different detection methods. (D) Classification of detections by confidence level:\n"
                "            High (≥7), Medium (4-6), and Low (1-3) confidence. Overall detection rate: {detection_rate:.1f}%.\n            "
            )
        }

        # Provide safe defaults to avoid KeyError when optional stats are missing
        defaults = {
            'detection_pairs': 0,
            'total_samples': 0,
            'total_pathogens': 0,
            'total_assessments': 0,
            'methods_analyzed': 0,
            'mean_pathogens_per_sample': 0.0,
            'detection_rate': 0.0,
        }
        merged = {**defaults, **stats}

        template = templates.get(figure_type)
        if not template:
            return "Figure caption not available."
        try:
            return template.format(**merged)
        except Exception:
            # As a last resort, return an unformatted generic caption
            return "Figure caption not available."
    
    def save_caption(self, figure_type, caption):
        """Save caption to file"""
        caption_file = os.path.join(self.output_dir, f"{figure_type}_caption.txt")
        with open(caption_file, 'w') as f:
            f.write(caption.strip())
    
    def generate_all_figures(self):
        """Generate all publication-ready figures"""
        print("Generating publication-ready figures...")
        
        figures = []
        
        # Generate all figure types
        if 'scores' in self.data:
            figures.append(self.create_multi_method_comparison())
            figures.append(self.create_sample_comparison())
            figures.append(self.create_summary_statistics())
        
        if 'detailed' in self.data:
            figures.append(self.create_method_performance_analysis())
        
        # Generate figure index
        self.create_figure_index(figures)
        
        print(f"All publication figures generated in: {self.output_dir}")
        return figures
    
    def create_figure_index(self, figures):
        """Create an index of all generated figures"""
        index_content = f"""
# PIGSTI Publication Figures Index
Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

## Available Figures

"""
        
        figure_info = [
            ("multi_method_comparison", "Multi-Method Pathogen Detection Analysis", "Comprehensive comparison of all detection methods"),
            ("method_performance_analysis", "Detection Method Performance Analysis", "Detailed analysis of individual method performance"),
            ("sample_comparison", "Sample-Level Pathogen Detection Analysis", "Sample-by-sample detection analysis"),
            ("summary_statistics", "PIGSTI Pipeline Summary Statistics", "Overall pipeline performance summary")
        ]
        
        for i, (fig_type, title, description) in enumerate(figure_info, 1):
            if any(fig_type in fig for fig in figures if fig):
                index_content += f"""
### Figure {i}: {title}
- **File**: `{fig_type}.png` / `{fig_type}.pdf`
- **Description**: {description}
- **Caption**: See `{fig_type}_caption.txt`

"""
        
        index_content += """
## File Formats
- **PNG**: High-resolution raster images (300 DPI) for presentations
- **PDF**: Vector graphics for publications and scalable use
- **Caption**: Text files with publication-ready captions

## Usage Notes
- All figures are generated at 300 DPI for publication quality
- Consistent color schemes and styling across all figures
- Captions are automatically generated and can be customized
- Figures are optimized for both print and digital publication
"""
        
        index_file = os.path.join(self.output_dir, "README.md")
        with open(index_file, 'w') as f:
            f.write(index_content)

def main():
    """Main function to generate publication figures"""
    generator = PublicationFigureGenerator()
    generator.generate_all_figures()

if __name__ == "__main__":
    main()
