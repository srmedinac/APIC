"""
APIC Feature Selection Pipeline
Part of: "An AI based Pathology Model to Predict Docetaxel Benefit in Prostate Cancer"

This script performs feature selection for the APIC classifier using elastic net 
penalized Cox regression on imaging features extracted from H&E-stained slides.
"""

import numpy as np
import pandas as pd
import warnings
from sklearn.feature_selection import VarianceThreshold
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import MinMaxScaler
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.util import Surv

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

class APICFeatureSelector:
    """Feature selection pipeline for APIC classifier development"""
    
    def __init__(self, l1_ratio=0.9, variance_threshold=0.01, correlation_threshold=0.85):
        """
        Initialize feature selector
        
        Args:
            l1_ratio: Elastic net mixing parameter (0=Ridge, 1=Lasso)
            variance_threshold: Minimum variance threshold for feature inclusion
            correlation_threshold: Maximum correlation threshold between features
        """
        self.l1_ratio = l1_ratio
        self.variance_threshold = variance_threshold
        self.correlation_threshold = correlation_threshold
        self.scaler = MinMaxScaler()
        self.imputer = SimpleImputer(strategy='median')
        
    def preprocess_features(self, df, feature_cols):
        """
        Preprocess imaging features: impute missing values and scale
        
        Args:
            df: Input dataframe
            feature_cols: List of feature column names
            
        Returns:
            Preprocessed feature matrix
        """
        # Handle infinite values
        feature_data = df[feature_cols].replace([np.inf, -np.inf], np.nan)
        
        # Impute missing values
        feature_data = pd.DataFrame(
            self.imputer.fit_transform(feature_data),
            columns=feature_cols,
            index=df.index
        )
        
        # Scale features to [0,1] range
        feature_data = pd.DataFrame(
            self.scaler.fit_transform(feature_data),
            columns=feature_cols,
            index=df.index
        )
        
        return feature_data
    
    def remove_low_variance_features(self, X):
        """Remove features with low variance"""
        selector = VarianceThreshold(threshold=self.variance_threshold)
        X_selected = selector.fit_transform(X)
        selected_features = X.columns[selector.get_support()]
        
        return pd.DataFrame(X_selected, columns=selected_features, index=X.index)
    
    def remove_correlated_features(self, X):
        """Remove highly correlated features"""
        corr_matrix = X.corr().abs()
        upper_triangle = corr_matrix.where(
            np.triu(np.ones(corr_matrix.shape), k=1).astype(bool)
        )
        
        to_drop = [col for col in upper_triangle.columns 
                  if any(upper_triangle[col] > self.correlation_threshold)]
        
        return X.drop(columns=to_drop)
    
    def select_features_bootstrap(self, X, y_surv, n_features=7, n_bootstrap=50, 
                                 alphas=None, random_state=42):
        """
        Select top features using bootstrap elastic net penalized Cox regression
        
        Args:
            X: Feature matrix
            y_surv: Survival data (structured array from sksurv)
            n_features: Number of top features to select
            n_bootstrap: Number of bootstrap iterations
            alphas: Regularization strengths to try
            random_state: Random seed for reproducibility
            
        Returns:
            List of selected feature names based on selection frequency
        """
        if alphas is None:
            alphas = [0.1, 0.01, 0.001, 0.0001]
        
        np.random.seed(random_state)
        feature_counts = pd.Series(0, index=X.columns)
        
        print(f"   - Running {n_bootstrap} bootstrap iterations...")
        
        for i in range(n_bootstrap):
            # Bootstrap sampling
            n_samples = len(X)
            bootstrap_idx = np.random.choice(n_samples, size=n_samples, replace=True)
            
            X_boot = X.iloc[bootstrap_idx]
            y_boot = y_surv[bootstrap_idx]
            
            # Fit elastic net Cox model
            cox_net = CoxnetSurvivalAnalysis(
                l1_ratio=self.l1_ratio,
                alphas=alphas,
                fit_baseline_model=True,
                normalize=False  # Already normalized
            )
            
            try:
                cox_net.fit(X_boot.values, y_boot)
                
                # Get non-zero coefficients from the most regularized model
                coefs = cox_net.coef_[0]
                selected_features = X_boot.columns[coefs != 0]
                
                # Count feature selections
                feature_counts[selected_features] += 1
                
            except Exception:
                # Skip failed iterations
                continue
            
            if (i + 1) % 10 == 0:
                print(f"     Completed {i + 1}/{n_bootstrap} iterations")
        
        # Select features based on selection frequency
        feature_frequencies = feature_counts / n_bootstrap
        top_features = feature_frequencies.nlargest(n_features).index.tolist()
        
        print(f"   - Feature selection frequencies:")
        for feature in top_features:
            freq = feature_frequencies[feature]
            print(f"     {feature:<30} selected in {freq:.1%} of bootstraps")
        
        return top_features, feature_frequencies
    
    def run_feature_selection(self, df, feature_cols, time_col='time', event_col='event', 
                            n_features=7, n_bootstrap=50, alphas=None, use_bootstrap=True):
        """
        Complete feature selection pipeline with bootstrap methodology
        
        Args:
            df: Input dataframe with features and outcomes
            feature_cols: List of feature column names
            time_col: Name of survival time column
            event_col: Name of event indicator column
            n_features: Number of features to select
            n_bootstrap: Number of bootstrap iterations
            alphas: Regularization parameters for elastic net
            use_bootstrap: Whether to use bootstrap feature selection
            
        Returns:
            Dictionary with selected features and preprocessing objects
        """
        print("=== APIC Feature Selection Pipeline ===\n")
        print(f"Input features: {len(feature_cols)}")
        print(f"Input samples: {len(df)}")
        print(f"Bootstrap iterations: {n_bootstrap}" if use_bootstrap else "Single model fit")
        
        # Step 1: Preprocess features
        print("\n1. Preprocessing features...")
        X_processed = self.preprocess_features(df, feature_cols)
        print(f"   - Imputed missing values and scaled to [0,1]")
        
        # Step 2: Remove low variance features
        print("2. Removing low variance features...")
        X_variance = self.remove_low_variance_features(X_processed)
        print(f"   - Removed {len(feature_cols) - len(X_variance.columns)} low variance features")
        print(f"   - Remaining features: {len(X_variance.columns)}")
        
        # Step 3: Remove correlated features
        print("3. Removing highly correlated features...")
        X_correlation = self.remove_correlated_features(X_variance)
        print(f"   - Removed {len(X_variance.columns) - len(X_correlation.columns)} correlated features")
        print(f"   - Remaining features: {len(X_correlation.columns)}")
        
        # Step 4: Create survival data
        y_surv = Surv.from_dataframe(event_col, time_col, df)
        
        # Step 5: Feature selection with bootstrap
        if use_bootstrap:
            print(f"4. Bootstrap feature selection (n={n_bootstrap})...")
            selected_features, feature_frequencies = self.select_features_bootstrap(
                X_correlation, y_surv, n_features, n_bootstrap, alphas
            )
        else:
            print(f"4. Single model feature selection...")
            selected_features, cox_model = self.select_features_elastic_net(
                X_correlation, y_surv, n_features, alphas
            )
            feature_frequencies = None
        
        print(f"\n=== Selected Features ===")
        for i, feature in enumerate(selected_features, 1):
            if use_bootstrap:
                freq = feature_frequencies[feature]
                print(f"{i:2d}. {feature:<30} (selected {freq:.1%} of bootstraps)")
            else:
                coef = cox_model.coef_[0][X_correlation.columns.get_loc(feature)]
                print(f"{i:2d}. {feature:<30} (coef: {coef:+.4f})")
        
        # Prepare results
        results = {
            'selected_features': selected_features,
            'feature_frequencies': feature_frequencies if use_bootstrap else None,
            'n_input_features': len(feature_cols),
            'n_after_variance': len(X_variance.columns),
            'n_after_correlation': len(X_correlation.columns),
            'n_selected': len(selected_features),
            'n_bootstrap': n_bootstrap if use_bootstrap else 1,
            'preprocessing': {
                'scaler': self.scaler,
                'imputer': self.imputer
            },
            'final_feature_matrix': X_correlation[selected_features]
        }
        
        print(f"\n=== Feature Selection Complete ===")
        if use_bootstrap:
            print(f"Selected {len(selected_features)} features using {n_bootstrap} bootstrap iterations")
        else:
            print(f"Selected {len(selected_features)} features using single model")
        
        return results