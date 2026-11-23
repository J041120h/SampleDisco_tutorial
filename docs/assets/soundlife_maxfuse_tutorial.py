# SampleDisco tutorial: unpaired scRNA-seq + spectral flow cytometry (Sound Life) with MaxFuse.
# Run from the folder that holds data/soundlife_RNA.h5ad and data/soundlife_FLOW.h5ad.
# Each "# %% Step" block is one section of the tutorial page.

# %% Step 1 - load the data
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

rna = ad.read_h5ad("data/soundlife_RNA.h5ad")    # raw counts
flow = ad.read_h5ad("data/soundlife_FLOW.h5ad")  # X = arcsinh(unmixed / 3000)
print(rna, flow, sep="\n")

# %% Step 2 - MaxFuse inputs
# Flow marker -> gene(s) encoding it. CD45, CD45RA and CD15 have no usable gene and stay flow-only.
LINKS = {
    "CD3": ["CD3E", "CD3D", "CD3G"], "CD19": ["CD19"], "CD56": ["NCAM1"], "CD14": ["CD14"],
    "HLA-DR": ["HLA-DRA", "HLA-DRB1"], "CD4": ["CD4"], "CD8": ["CD8A", "CD8B"], "CD25": ["IL2RA"],
    "CD127": ["IL7R"], "CCR7": ["CCR7"], "TCRab": ["TRAC"], "IgD": ["IGHD"], "CD27": ["CD27"],
    "CD38": ["CD38"], "CD16": ["FCGR3A"], "CD11b": ["ITGAM"], "CD11c": ["ITGAX"], "CD123": ["IL3RA"],
    "CD304": ["NRP1"], "CD141": ["THBD"], "CD34": ["CD34"],
}
pairs = [(g, p) for p, genes in LINKS.items() for g in genes
         if g in rna.var_names and p in flow.var_names]
rna_link = rna[:, [g for g, _ in pairs]].X.toarray()
flow_link = flow[:, [p for _, p in pairs]].X
keep = (rna_link.std(0) > 1e-5) & (flow_link.std(0) > 1e-5)   # drop constant pairs
pairs = [pr for pr, k in zip(pairs, keep) if k]
print(f"{len(pairs)} linked gene-protein pairs")

def zscore(x):
    return ((x - x.mean(0)) / x.std(0)).astype(np.float32)

# shared features: linked genes normalized to the flow scale, log1p, z-scored
shared = ad.AnnData(rna[:, [g for g, _ in pairs]].X.astype(np.float32))
target = (np.median(np.asarray(shared.X.sum(1))) + np.median(flow_link[:, keep].sum(1))) / 2
sc.pp.normalize_total(shared, target_sum=target)
sc.pp.log1p(shared)
rna_shared = zscore(shared.X.toarray())
flow_z = zscore(flow.X)
flow_shared = flow_z[:, [flow.var_names.get_loc(p) for _, p in pairs]]

# active features: RNA highly variable genes, all flow markers
act = rna.copy()
sc.pp.normalize_total(act)
sc.pp.log1p(act)
sc.pp.highly_variable_genes(act)
rna_active = zscore(act[:, act.var["highly_variable"]].X.toarray())
flow_active = flow_z
print(f"shared {rna_shared.shape[1]} features, active RNA {rna_active.shape[1]} HVGs / flow {flow_active.shape[1]} markers")

# %% Step 3 - run MaxFuse
import maxfuse as mf

n_shared = rna_shared.shape[1] - 1
fusor = mf.model.Fusor(shared_arr1=rna_shared, shared_arr2=flow_shared,
                       active_arr1=rna_active, active_arr2=flow_active)
fusor.split_into_batches(max_outward_size=5000, matching_ratio=3, metacell_size=2, seed=42)
fusor.construct_graphs(n_neighbors1=15, n_neighbors2=15, svd_components1=30, svd_components2=None,
                       resolution1=2, resolution2=2, resolution_tol=0.1, leiden_seed=42)
fusor.find_initial_pivots(wt1=0.7, wt2=0.7, svd_components1=min(25, n_shared),
                          svd_components2=min(20, n_shared))
fusor.refine_pivots(wt1=0.7, wt2=0.7, svd_components1=30, svd_components2=None, cca_components=20,
                    n_iters=3, randomized_svd=False, svd_runs=1)
fusor.filter_bad_matches(target="pivot", filter_prop=0.3)
rna_cca, flow_cca = fusor.get_embedding(active_arr1=fusor.active_arr1, active_arr2=fusor.active_arr2)
print("joint embedding:", rna_cca.shape, flow_cca.shape)

# %% Step 4 - flow cell types and the RNA + flow union
from scipy import sparse
from sklearn.neighbors import NearestNeighbors
from sampledisco.utils.harmony_compat import harmonize_embedding
from sampledisco.preparation.multi_omics_batch_correction import harmonize_xglue

# each flow event takes the majority RNA cell type among its 15 nearest RNA cells in the joint space
_, nbr = NearestNeighbors(n_neighbors=15).fit(rna_cca).kneighbors(flow_cca)
codes, names = pd.factorize(rna.obs["cell_type"].astype(str))
flow_ct = names[[np.bincount(row, minlength=len(names)).argmax() for row in codes[nbr]]]

# sanity check: transferred lineage vs the flow marker gates
to_lineage = {"T cell": "T", "B cell": "B", "NK cell": "NK", "Monocyte": "Mono", "DC": "DC"}
l1 = rna.obs.groupby("cell_type", observed=True)["cell_type_l1"].first().astype(str)
flow_lin = pd.Series(flow_ct).map(l1).map(to_lineage).fillna("other").to_numpy()
print(f"lineage agreement with flow gates: {(flow_lin == flow.obs['lineage'].to_numpy()).mean():.3f}")

meta = ["sample", "batch", "subject", "visit", "age_group", "age_at_draw", "sex", "cmv"]
obs = pd.concat([
    rna.obs[meta].astype({c: str for c in meta if c != "age_at_draw"})
       .assign(modality="RNA", cell_type=rna.obs["cell_type"].astype(str).to_numpy()),
    flow.obs[meta].astype({c: str for c in meta if c != "age_at_draw"})
        .assign(modality="FLOW", cell_type=np.asarray(flow_ct)),
])
for c in obs.columns.drop("age_at_draw"):
    obs[c] = obs[c].astype("category")
union = ad.AnnData(X=sparse.csr_matrix((len(obs), 0), dtype=np.float32), obs=obs)
union.obsm["X_maxfuse"] = np.vstack([rna_cca, flow_cca]).astype(np.float32)

# Z_rmd: technical batch removed, sample differences kept (displacement view)
union.obsm["Z_rmd"] = harmonize_embedding(union.obsm["X_maxfuse"], union.obs, batch_key=["batch"],
                                          max_iter_harmony=50, use_gpu=False, seed=42)
# Z_comp: batch and sample removed (composition / clustering view)
harmonize_xglue(union, sample_col="sample", batch_col="batch", in_key="Z_rmd", out_key="Z_comp",
                use_gpu=False, max_iter=50, random_state=42)
print(union)

# %% Step 5 - sample embedding
from sampledisco.sample_embedding import compute_sample_embedding

union = compute_sample_embedding(
    union,
    output_dir="sampledisco_soundlife_output/multiomics",
    use_gpu=False,
    sample_col="sample",
    celltype_col="cell_type",
    comp_emb_key="Z_comp",
    rmd_emb_key="Z_rmd",
    modality_col="modality",
    batch_col="batch",
    medium_K=120,
    fine_K=300,
    rmd_dim_per_cluster=8,
    rmd_weight=0.60,
    pca_components=10,
    batch_method="harmony",
    save=True,
)
print(union.uns["X_DR_sample"].iloc[:4, :4])

# %% Step 6 - what the sample embedding captures
from sampledisco.sample_embedding.sample_embedding import build_sample_adata
from sampledisco.sample_association.association import run_dimension_association_analysis
from sampledisco.visualization.multi_omics_visualization import visualize_multimodal_embedding

sample_adata = build_sample_adata(union, sample_col="sample", modality_col="modality")
assoc = run_dimension_association_analysis(
    pseudo_adata=sample_adata,
    output_dir="sampledisco_soundlife_output/multiomics/sample_association",
    continuous_cols=["age_at_draw"],
    categorical_cols=["age_group", "sex", "cmv"],
    n_permutations=999,
    sample_col="sample",
    random_state=42,
)
visualize_multimodal_embedding(
    adata=sample_adata,
    visualization_grouping_column=["age_group", "modality"],
    output_dir="sampledisco_soundlife_output/multiomics",   # plots go to <output_dir>/visualization/
    show_default=True,
)
