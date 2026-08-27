"""
Script to evaluate alternative E-score functions on KrakenUniq data
Compares discrimination ability between true and false positives
"""

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from alternative_escore_functions import (
    e_score_dexp, e_score_guellil, e_score_sigmoid, e_score_logistic,
    e_score_adaptive_exp, e_score_with_dup_penalty, e_score_ratio_focused,
    e_score_balanced, e_score_tanh, e_score_piecewise, e_score_robust,
    e_score_multiplicative_penalty, compare_scores, get_recommended_functions
)

def load_kraken_data(escore_file):
    """Load KrakenUniq Escore CSV file"""
    df = pd.read_csv(escore_file)
    df.columns = df.columns.str.strip()
    
    # Ensure required columns exist
    required = ['uniq_kmers', 'reads', 'cov']
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    
    # Handle duplication rate if available
    if 'dup' in df.columns:
        df['dup_rate'] = df['dup'] / (df['reads'] + 1e-10)
    else:
        df['dup_rate'] = 0.0
    
    return df


def simulate_true_false_positives(df, n_samples=100):
    """
    Simulate true and false positive scenarios based on Maxime Borry's blog post
    This helps evaluate discrimination without needing labeled data
    """
    results = []
    
    # Scenario A: High duplication, low coverage (FALSE POSITIVE)
    # Many reads are duplicated, probably from conserved regions
    for _ in range(n_samples):
        reads = np.random.randint(10, 300)
        kmers = int(reads / np.random.randint(2, 10))
        cov = kmers / 1000000  # Assuming ~1M genome size
        dup_rate = np.random.uniform(0.5, 0.9)
        results.append({
            'uniq_kmers': kmers,
            'reads': reads,
            'cov': cov,
            'dup_rate': dup_rate,
            'label': 'false_positive',
            'scenario': 'A'
        })
    
    # Scenario B: Low duplication, low coverage (TRUE POSITIVE)
    # Very few reads duplicated, more unique k-mers than reads
    for _ in range(n_samples):
        reads = np.random.randint(10, 300)
        kmers = min(1000000 - 35 + 1, int(reads * np.random.randint(2, 11)))
        cov = kmers / 1000000
        dup_rate = np.random.uniform(0.0, 0.3)
        results.append({
            'uniq_kmers': kmers,
            'reads': reads,
            'cov': cov,
            'dup_rate': dup_rate,
            'label': 'true_positive',
            'scenario': 'B'
        })
    
    # Scenario C: Low duplication, higher coverage (TRUE POSITIVE)
    for _ in range(n_samples):
        reads = np.random.randint(10000, 100000)
        kmers = min(1000000 - 35 + 1, int(reads * np.random.randint(2, 11)))
        cov = kmers / 1000000
        dup_rate = np.random.uniform(0.0, 0.3)
        results.append({
            'uniq_kmers': kmers,
            'reads': reads,
            'cov': cov,
            'dup_rate': dup_rate,
            'label': 'true_positive',
            'scenario': 'C'
        })
    
    return pd.DataFrame(results)


def evaluate_functions(df, functions_dict):
    """Evaluate scoring functions and return discrimination metrics"""
    results_df = df.copy()
    
    # Calculate scores
    for name, func in functions_dict.items():
        try:
            if 'dup' in str(func.__code__.co_varnames) and 'dup_rate' in df.columns:
                results_df[f'score_{name}'] = df.apply(
                    lambda row: func(row['uniq_kmers'], row['reads'], row['cov'], row['dup_rate']),
                    axis=1
                )
            else:
                results_df[f'score_{name}'] = df.apply(
                    lambda row: func(row['uniq_kmers'], row['reads'], row['cov']),
                    axis=1
                )
        except Exception as e:
            print(f"Warning: Could not calculate {name}: {e}")
            results_df[f'score_{name}'] = np.nan
    
    # Calculate discrimination metrics if labels available
    if 'label' in df.columns:
        metrics = {}
        tp_mask = df['label'] == 'true_positive'
        fp_mask = df['label'] == 'false_positive'
        
        for name in functions_dict.keys():
            score_col = f'score_{name}'
            if score_col not in results_df.columns:
                continue
            
            tp_scores = results_df.loc[tp_mask, score_col].dropna()
            fp_scores = results_df.loc[fp_mask, score_col].dropna()
            
            if len(tp_scores) > 0 and len(fp_scores) > 0:
                mean_diff = tp_scores.mean() - fp_scores.mean()
                std_pooled = np.sqrt(tp_scores.var() + fp_scores.var())
                separation = mean_diff / (std_pooled + 1e-10) if std_pooled > 0 else 0
                
                # Calculate approximate AUC (area under ROC curve)
                threshold = fp_scores.median()
                tp_above_threshold = (tp_scores > threshold).sum()
                auc_approx = tp_above_threshold / len(tp_scores) if len(tp_scores) > 0 else 0
                
                metrics[name] = {
                    'mean_tp': tp_scores.mean(),
                    'mean_fp': fp_scores.mean(),
                    'mean_diff': mean_diff,
                    'separation': separation,
                    'auc_approx': auc_approx,
                    'median_tp': tp_scores.median(),
                    'median_fp': fp_scores.median(),
                }
        
        return results_df, pd.DataFrame(metrics).T
    else:
        return results_df, None


def plot_comparison(results_df, metrics_df, output_file=None):
    """Create comparison plots"""
    if metrics_df is None:
        print("No metrics available for plotting (missing labels)")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    
    # 1. Mean difference between TP and FP
    ax1 = axes[0, 0]
    metrics_df_sorted = metrics_df.sort_values('mean_diff', ascending=False)
    ax1.barh(range(len(metrics_df_sorted)), metrics_df_sorted['mean_diff'])
    ax1.set_yticks(range(len(metrics_df_sorted)))
    ax1.set_yticklabels(metrics_df_sorted.index)
    ax1.set_xlabel('Mean Score Difference (TP - FP)')
    ax1.set_title('Discrimination: Mean Difference')
    ax1.axvline(x=0, color='red', linestyle='--', alpha=0.5)
    ax1.grid(axis='x', alpha=0.3)
    
    # 2. Separation metric
    ax2 = axes[0, 1]
    metrics_df_sorted2 = metrics_df.sort_values('separation', ascending=False)
    ax2.barh(range(len(metrics_df_sorted2)), metrics_df_sorted2['separation'])
    ax2.set_yticks(range(len(metrics_df_sorted2)))
    ax2.set_yticklabels(metrics_df_sorted2.index)
    ax2.set_xlabel('Separation (mean_diff / pooled_std)')
    ax2.set_title('Discrimination: Separation Score')
    ax2.grid(axis='x', alpha=0.3)
    
    # 3. Approximate AUC
    ax3 = axes[1, 0]
    metrics_df_sorted3 = metrics_df.sort_values('auc_approx', ascending=False)
    ax3.barh(range(len(metrics_df_sorted3)), metrics_df_sorted3['auc_approx'])
    ax3.set_yticks(range(len(metrics_df_sorted3)))
    ax3.set_yticklabels(metrics_df_sorted3.index)
    ax3.set_xlabel('Approximate AUC')
    ax3.set_title('Discrimination: Approximate AUC')
    ax3.axvline(x=0.5, color='red', linestyle='--', alpha=0.5, label='Random (0.5)')
    ax3.legend()
    ax3.grid(axis='x', alpha=0.3)
    
    # 4. Score distributions (top 3 functions)
    ax4 = axes[1, 1]
    top_3 = metrics_df.nlargest(3, 'separation').index
    if len(top_3) > 0:
        for name in top_3:
            score_col = f'score_{name}'
            if score_col in results_df.columns:
                tp_scores = results_df[results_df['label'] == 'true_positive'][score_col].dropna()
                fp_scores = results_df[results_df['label'] == 'false_positive'][score_col].dropna()
                ax4.hist(tp_scores, alpha=0.5, label=f'{name} (TP)', bins=30, density=True)
                ax4.hist(fp_scores, alpha=0.5, label=f'{name} (FP)', bins=30, density=True)
        ax4.set_xlabel('Score')
        ax4.set_ylabel('Density')
        ax4.set_title('Score Distributions (Top 3 Functions)')
        ax4.legend()
        ax4.grid(alpha=0.3)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    else:
        plt.show()


def main():
    """Main evaluation function"""
    if len(sys.argv) < 2:
        print("Usage: python evaluate_escore_functions.py <escore_file.csv> [output_dir]")
        print("\nIf no escore file provided, will use simulated data")
        sys.exit(1)
    
    escore_file = sys.argv[1] if len(sys.argv) > 1 else None
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("results/escore_evaluation")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load or simulate data
    if escore_file and Path(escore_file).exists():
        print(f"Loading data from {escore_file}")
        df = load_kraken_data(escore_file)
        # If no labels, add simulated scenarios for comparison
        if 'label' not in df.columns:
            print("No labels found, adding simulated true/false positive scenarios")
            sim_df = simulate_true_false_positives(df.head(100), n_samples=50)
            df = pd.concat([df, sim_df], ignore_index=True)
    else:
        print("No valid escore file, using simulated data")
        df = simulate_true_false_positives(None, n_samples=100)
    
    # Get recommended functions
    functions = get_recommended_functions()
    
    # Add duplication-based functions if dup_rate available
    if 'dup_rate' in df.columns:
        functions['with_dup_penalty'] = lambda k, r, c, d: e_score_with_dup_penalty(k, r, c, d)
        functions['multiplicative'] = lambda k, r, c, d: e_score_multiplicative_penalty(k, r, c, d, min_reads=10)
    
    print(f"\nEvaluating {len(functions)} scoring functions...")
    
    # Evaluate
    results_df, metrics_df = evaluate_functions(df, functions)
    
    # Save results
    results_df.to_csv(output_dir / "evaluation_results.csv", index=False)
    if metrics_df is not None:
        metrics_df.to_csv(output_dir / "discrimination_metrics.csv")
        print("\nDiscrimination Metrics:")
        print(metrics_df.sort_values('separation', ascending=False))
        
        # Create plots
        plot_comparison(results_df, metrics_df, output_dir / "escore_comparison.png")
    
    # Recommendations
    if metrics_df is not None:
        print("\n" + "="*60)
        print("RECOMMENDATIONS:")
        print("="*60)
        best_separation = metrics_df.nlargest(1, 'separation').index[0]
        best_auc = metrics_df.nlargest(1, 'auc_approx').index[0]
        best_diff = metrics_df.nlargest(1, 'mean_diff').index[0]
        
        print(f"\nBest Separation Score: {best_separation}")
        print(f"  - Separation: {metrics_df.loc[best_separation, 'separation']:.3f}")
        print(f"  - Mean Difference: {metrics_df.loc[best_separation, 'mean_diff']:.3f}")
        
        print(f"\nBest AUC Approximation: {best_auc}")
        print(f"  - AUC: {metrics_df.loc[best_auc, 'auc_approx']:.3f}")
        
        print(f"\nBest Mean Difference: {best_diff}")
        print(f"  - Mean Difference: {metrics_df.loc[best_diff, 'mean_diff']:.3f}")
        
        print("\nTop 3 Functions Overall (by separation):")
        top_3 = metrics_df.nlargest(3, 'separation')
        for i, (name, row) in enumerate(top_3.iterrows(), 1):
            print(f"  {i}. {name}: separation={row['separation']:.3f}, "
                  f"auc={row['auc_approx']:.3f}, diff={row['mean_diff']:.3f}")


if __name__ == "__main__":
    main()

