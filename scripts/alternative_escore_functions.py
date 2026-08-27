"""
Alternative E-score functions for improved true/false positive discrimination
Based on Maxime Borry's blog post: https://maximeborry.com/post/kraken-uniq/

This module provides multiple scoring functions that can be tested and compared
for their ability to discriminate true positives from false positives in KrakenUniq outputs.
"""

import numpy as np
import pandas as pd


# ============================================================================
# EXISTING FUNCTIONS (for comparison)
# ============================================================================

def double_exp(x, a=1.3, b=18):
    """Double exponential function used in current Escore"""
    try:
        return a ** (b * x)
    except OverflowError:
        return float("inf")


def e_score_guellil(nb_kmer, nb_read, cov):
    """Original Guellil et al. E-value: (K/R) * C"""
    if nb_read == 0:
        return 0
    return (nb_kmer / nb_read) * cov


def e_score_dexp(nb_kmer, nb_read, cov):
    """Current double-exponential E-score: (K/R) * dexp(C)"""
    if nb_read == 0:
        return 0
    return (nb_kmer / nb_read) * double_exp(cov)


# ============================================================================
# ALTERNATIVE SCORING FUNCTIONS
# ============================================================================

def sigmoid(x, k=10, x0=0.5):
    """
    Sigmoid function: S(x) = 1 / (1 + exp(-k*(x-x0)))
    Provides smooth transition, less aggressive than exponential
    """
    try:
        return 1 / (1 + np.exp(-k * (x - x0)))
    except OverflowError:
        return 1.0 if x > x0 else 0.0


def e_score_sigmoid(nb_kmer, nb_read, cov, k=10, x0=0.5):
    """
    E-score using sigmoid function instead of double exponential
    Parameters:
        k: steepness of sigmoid (default 10)
        x0: midpoint of sigmoid (default 0.5)
    """
    if nb_read == 0:
        return 0
    return (nb_kmer / nb_read) * sigmoid(cov, k=k, x0=x0)


def e_score_logistic(nb_kmer, nb_read, cov, L=100, k=20, x0=0.3):
    """
    Logistic function: L / (1 + exp(-k*(x-x0)))
    Similar to sigmoid but with adjustable maximum value
    """
    if nb_read == 0:
        return 0
    try:
        logistic = L / (1 + np.exp(-k * (cov - x0)))
    except OverflowError:
        logistic = L if cov > x0 else 0
    return (nb_kmer / nb_read) * logistic


def e_score_adaptive_exp(nb_kmer, nb_read, cov, base=1.3, power_mult=18, min_cov=0.01):
    """
    Adaptive exponential: uses exponential but with floor to prevent
    extreme values at very low coverage
    """
    if nb_read == 0:
        return 0
    # Apply exponential only above minimum coverage threshold
    if cov < min_cov:
        # Linear scaling for very low coverage
        exp_term = base ** (power_mult * min_cov) * (cov / min_cov)
    else:
        exp_term = base ** (power_mult * cov)
    return (nb_kmer / nb_read) * exp_term


def e_score_with_dup_penalty(nb_kmer, nb_read, cov, dup_rate):
    """
    E-score that penalizes high duplication rates
    False positives often have high duplication (read stacking)
    """
    if nb_read == 0:
        return 0
    # Penalty factor: (1 - dup_rate) means lower score for higher duplication
    # dup_rate typically ranges from 0 to 1
    penalty = 1.0 - min(dup_rate, 0.95)  # Cap penalty to avoid division issues
    return (nb_kmer / nb_read) * double_exp(cov) * penalty


def e_score_ratio_focused(nb_kmer, nb_read, cov, alpha=2.0):
    """
    Emphasizes the k-mer to read ratio more strongly
    True positives typically have higher K/R ratios
    """
    if nb_read == 0:
        return 0
    kmer_ratio = nb_kmer / nb_read
    # Apply power to ratio to emphasize it
    return (kmer_ratio ** alpha) * double_exp(cov)


def e_score_balanced(nb_kmer, nb_read, cov, w1=0.6, w2=0.4):
    """
    Balanced score combining ratio and coverage components
    Allows tuning the relative importance of each component
    """
    if nb_read == 0:
        return 0
    kmer_ratio = nb_kmer / nb_read
    # Weighted combination
    return w1 * kmer_ratio * double_exp(cov) + w2 * kmer_ratio * cov


def e_score_tanh(nb_kmer, nb_read, cov, scale=2.0):
    """
    Hyperbolic tangent function: tanh(scale * x)
    Provides smooth saturation, less extreme than exponential
    """
    if nb_read == 0:
        return 0
    return (nb_kmer / nb_read) * np.tanh(scale * cov)


def e_score_piecewise(nb_kmer, nb_read, cov, threshold=0.3):
    """
    Piecewise function: different behavior below and above coverage threshold
    Can help distinguish low-coverage true positives from false positives
    """
    if nb_read == 0:
        return 0
    kmer_ratio = nb_kmer / nb_read
    
    if cov < threshold:
        # For low coverage: use linear scaling (less aggressive)
        return kmer_ratio * (cov / threshold) * double_exp(threshold)
    else:
        # For higher coverage: use full exponential
        return kmer_ratio * double_exp(cov)


def e_score_multiplicative_penalty(nb_kmer, nb_read, cov, dup_rate, min_reads=10):
    """
    Multi-factor score combining multiple signals
    Penalizes: high duplication, low read counts
    Rewards: high k-mer diversity, good coverage
    """
    if nb_read == 0:
        return 0
    
    kmer_ratio = nb_kmer / nb_read
    
    # Read count penalty (very low read counts are suspicious)
    read_penalty = min(1.0, nb_read / min_reads)
    
    # Duplication penalty
    dup_penalty = 1.0 - min(dup_rate, 0.9)
    
    # Coverage boost
    cov_boost = double_exp(cov)
    
    return kmer_ratio * cov_boost * read_penalty * dup_penalty


def e_score_robust(nb_kmer, nb_read, cov, dup_rate=None):
    """
    Robust score designed to handle edge cases better
    Uses log transformation to reduce extreme values
    """
    if nb_read == 0 or nb_kmer == 0:
        return 0
    
    kmer_ratio = nb_kmer / nb_read
    
    # Log transform of coverage to reduce extreme values
    # Add small epsilon to avoid log(0)
    log_cov = np.log1p(cov * 100)  # Scale coverage for better log behavior
    
    base_score = kmer_ratio * log_cov
    
    # Apply duplication penalty if available
    if dup_rate is not None:
        base_score *= (1.0 - min(dup_rate, 0.9))
    
    return base_score


# ============================================================================
# COMPARISON AND EVALUATION FUNCTIONS
# ============================================================================

def compare_scores(df, score_functions, true_positive_label='true_positive'):
    """
    Compare multiple scoring functions on a dataframe
    
    Parameters:
        df: DataFrame with columns: uniq_kmers, reads, cov, (optional: dup, label)
        score_functions: dict of {name: function} where function takes (nb_kmer, nb_read, cov, ...)
        true_positive_label: column name or value indicating true positives
    
    Returns:
        DataFrame with scores for each function and discrimination metrics
    """
    results = df.copy()
    
    # Calculate scores for each function
    for name, func in score_functions.items():
        if 'dup' in df.columns and 'dup_rate' in str(func.__code__.co_varnames):
            # Function uses duplication rate
            results[f'score_{name}'] = df.apply(
                lambda row: func(row['uniq_kmers'], row['reads'], row['cov'], 
                                row.get('dup', 0)), axis=1
            )
        else:
            # Function doesn't use duplication
            results[f'score_{name}'] = df.apply(
                lambda row: func(row['uniq_kmers'], row['reads'], row['cov']), axis=1
            )
    
    # Calculate discrimination metrics if labels are available
    if true_positive_label in df.columns or isinstance(true_positive_label, str):
        if true_positive_label in df.columns:
            tp_mask = df[true_positive_label] == True
        else:
            tp_mask = df['label'] == true_positive_label
        
        fp_mask = ~tp_mask
        
        discrimination_metrics = {}
        for name in score_functions.keys():
            score_col = f'score_{name}'
            tp_scores = results.loc[tp_mask, score_col]
            fp_scores = results.loc[fp_mask, score_col]
            
            if len(tp_scores) > 0 and len(fp_scores) > 0:
                discrimination_metrics[name] = {
                    'mean_tp': tp_scores.mean(),
                    'mean_fp': fp_scores.mean(),
                    'mean_diff': tp_scores.mean() - fp_scores.mean(),
                    'median_tp': tp_scores.median(),
                    'median_fp': fp_scores.median(),
                    'separation': (tp_scores.mean() - fp_scores.mean()) / 
                                (tp_scores.std() + fp_scores.std() + 1e-10),
                    'auc_approx': len(tp_scores[tp_scores > fp_scores.median()]) / len(tp_scores)
                }
        
        return results, pd.DataFrame(discrimination_metrics).T
    
    return results, None


# ============================================================================
# RECOMMENDED FUNCTION FOR TESTING
# ============================================================================

def get_recommended_functions():
    """
    Returns a dictionary of recommended scoring functions to test
    """
    return {
        'current_dexp': lambda k, r, c: e_score_dexp(k, r, c),
        'guellil_original': lambda k, r, c: e_score_guellil(k, r, c),
        'sigmoid_k10': lambda k, r, c: e_score_sigmoid(k, r, c, k=10, x0=0.5),
        'sigmoid_k20': lambda k, r, c: e_score_sigmoid(k, r, c, k=20, x0=0.3),
        'logistic': lambda k, r, c: e_score_logistic(k, r, c, L=100, k=20, x0=0.3),
        'adaptive_exp': lambda k, r, c: e_score_adaptive_exp(k, r, c, min_cov=0.01),
        'ratio_focused': lambda k, r, c: e_score_ratio_focused(k, r, c, alpha=2.0),
        'balanced': lambda k, r, c: e_score_balanced(k, r, c, w1=0.6, w2=0.4),
        'tanh': lambda k, r, c: e_score_tanh(k, r, c, scale=2.0),
        'piecewise': lambda k, r, c: e_score_piecewise(k, r, c, threshold=0.3),
        'robust': lambda k, r, c: e_score_robust(k, r, c),
    }


if __name__ == "__main__":
    # Example usage
    print("Alternative E-score functions for KrakenUniq discrimination")
    print("=" * 60)
    print("\nAvailable functions:")
    print("1. e_score_sigmoid - Smooth sigmoid transition")
    print("2. e_score_logistic - Adjustable logistic function")
    print("3. e_score_adaptive_exp - Exponential with floor")
    print("4. e_score_with_dup_penalty - Penalizes duplication")
    print("5. e_score_ratio_focused - Emphasizes K/R ratio")
    print("6. e_score_balanced - Weighted combination")
    print("7. e_score_tanh - Hyperbolic tangent")
    print("8. e_score_piecewise - Different behavior by coverage")
    print("9. e_score_multiplicative_penalty - Multi-factor")
    print("10. e_score_robust - Log-transformed robust score")
    print("\nUse compare_scores() to evaluate on your data!")

