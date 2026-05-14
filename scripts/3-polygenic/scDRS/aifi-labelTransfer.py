import scanpy as sc
import numpy as np
import scipy.sparse
import celltypist
import time

def preprocess_for_celltypist(adata, dataset_name="Dataset"):
    print(f"--- Processing {dataset_name} ---")

    # Helper to check if matrix contains integers (raw counts)
    def is_integer_matrix(matrix):
        if scipy.sparse.issparse(matrix):
            # Check a sample of the data to save time
            values = matrix.data[:100] 
        else:
            values = matrix.flat[:100]
        # Allow a tiny tolerance for floating point errors in storage
        return np.all(np.abs(values - np.round(values)) < 1e-6)

    # --- STEP 1: Ensure we have Raw Counts ---
    if is_integer_matrix(adata.X):
        print(f"✓ {dataset_name} detected as raw integer counts.")
    
    else:
        print(f"! {dataset_name} .X is not raw (contains floats). Attempting to retrieve raw counts...")
        
        # Priority A: Check for a 'counts' layer (Standard Best Practice)
        if 'counts' in adata.layers:
            print(f"  -> Found 'counts' layer. Restoring...")
            adata.X = adata.layers['counts'].copy()
            
        # Priority B: Check adata.raw (Legacy/Alternative Practice)
        elif adata.raw is not None:
            print(f"  -> Found adata.raw. Reverting using raw.to_adata()...")
            # raw.to_adata() creates a new object with the raw matrix as .X
            adata = adata.raw.to_adata()
        
        # Check again if we successfully retrieved integers
        if is_integer_matrix(adata.X):
             print(f"✓ Successfully restored raw counts for {dataset_name}.")
        else:
            print(f"⚠️ WARNING: Could not find integer counts in .layers['counts'] or .raw.")
            print(f"   Proceeding with current data, but results may be invalid if data was already log-transformed.")

    # --- STEP 2: Normalize and Log Transform ---
    
    # 1. Save raw counts to a layer if not already there (safety net)
    if 'counts' not in adata.layers:
        adata.layers['counts'] = adata.X.copy()

    # 2. Normalize to 10,000 (CP10k)
    print(f"  -> Normalizing to target sum 10,000...")
    sc.pp.normalize_total(adata, target_sum=1e4)

    # 3. Log-transform (log1p)
    print(f"  -> Log-transforming...")
    sc.pp.log1p(adata)
    
    print(f"✓ {dataset_name} is ready for Label Transfer.\n")
    return adata

# --- Usage ---

# 1. Load the Query Data
print("Loading query dataset...")
adata_query = sc.read_h5ad('data/immune_health_atlas.h5ad')

# 2. Preprocess the Query (CRITICAL STEP)
# This ensures your query has the same normalization (log1p CP10k) as the training data
adata_query = preprocess_for_celltypist(adata_query, "Query")

# 3. Run Prediction using the saved model
# You can pass the path to the .pkl file directly
print("Loading model and predicting labels...")
model_path = 'output/wg2_scpred_model.pkl'

predictions = celltypist.annotate(
    adata_query, 
    model=model_path, 
    majority_voting=True  # Recommended: Refines labels based on local cell clusters
)

# 4. Integrate results back into AnnData
adata_query = predictions.to_adata()

print("Prediction complete.")
print(adata_query.obs[['predicted_labels', 'over_clustering', 'majority_voting', 'conf_score']].head())

# 5. Save the annotated query
adata_query.write_h5ad('output/immune_health_atlas_annotated.h5ad')
