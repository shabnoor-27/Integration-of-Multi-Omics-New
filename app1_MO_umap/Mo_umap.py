import os
import tempfile
import requests
import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.graph_objects as go
from pyvis.network import Network
from gseapy import enrichr
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import umap


# -----------------------------
# App Configuration
# -----------------------------

st.image("https://raw.githubusercontent.com/priyadarshinikp1/Multiomics-Integrator-app/main/logo.png", width=200)
st.title("Multi-Omics Statistical Integration with UMAP Vizzhy App")

with st.sidebar:
    st.markdown("---")
    st.markdown("**Developed by: PRIYADARSHINI**")
    st.markdown("[LinkedIn](https://www.linkedin.com/in/priyadarshini24) | [GitHub](https://github.com/priyadarshinikp1)")



# -----------------------------
# Sidebar Configuration
# -----------------------------
st.sidebar.header("⚙️ Filter Settings")
cadd_thresh = float(st.sidebar.text_input("Min CADD Score (Genomics)", value="20"))
logfc_thresh = float(st.sidebar.text_input("Min |logFC| (Transcriptomics)", value="1"))
t_pval_thresh = float(st.sidebar.text_input("Max p-value (Transcriptomics)", value="0.05"))
p_intensity_thresh = float(st.sidebar.text_input("Min Intensity (Proteomics)", value="1000"))
n_clusters = st.sidebar.slider("Number of Clusters", min_value=2, max_value=10, value=4)

# -----------------------------
# File Upload Section
# -----------------------------
st.header("📁 Upload Omics Data")
genomics = st.file_uploader("Upload Genomics CSV", type="csv")
transcriptomics = st.file_uploader("Upload Transcriptomics CSV", type="csv")
proteomics = st.file_uploader("Upload Proteomics CSV", type="csv")

# -----------------------------
# Helper Functions
# -----------------------------
def extract_uniprot_ids(protein_series):
    ids = set()
    for entry in protein_series.dropna():
        for pid in str(entry).split(";"):
            if pid.strip():
                ids.add(pid.strip())
    return ids

def map_uniprot_to_gene(uniprot_ids):
    mapping = {}
    ids = list(uniprot_ids)
    for i in range(0, len(ids), 100):
        chunk = ids[i:i+100]
        query = " OR ".join([f"accession:{id_}" for id_ in chunk])
        url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&fields=accession,gene_names&format=tsv"
        try:
            r = requests.get(url)
            if r.status_code == 200:
                lines = r.text.strip().split('\n')[1:]
                for line in lines:
                    acc, genes = line.split('\t')
                    mapping[acc] = genes.split()[0] if genes else acc
        except:
            pass
    return mapping

# -----------------------------
# Main Processing
# -----------------------------
if genomics and transcriptomics and proteomics:
    gdf = pd.read_csv(genomics)
    tdf = pd.read_csv(transcriptomics)
    pdf = pd.read_csv(proteomics)

    gdf_filtered = gdf[gdf['CADD'] >= cadd_thresh]
    tdf_filtered = tdf[(tdf['p_value'] <= t_pval_thresh) ]
    # Filter Transcriptomics data based on logFC threshold
    if logfc_thresh > 0:
    # For positive logFC threshold, filter for values >= logFC_thresh
         tdf_filtered = tdf[tdf['logFC'] >= logfc_thresh]
    elif logfc_thresh < 0:
    # For negative logFC threshold, filter for values < logFC_thresh
         tdf_filtered = tdf[tdf['logFC'] < logfc_thresh]
    else:
    # If logFC_thresh is 0, no filtering needed
         tdf_filtered = tdf

    pdf_filtered = pdf[pdf['Intensity'] >= p_intensity_thresh]

    st.subheader("🔍 Filtered Data Preview")
    top_n = st.slider("Show Top N Rows", 5, 50, 10)
    st.write("**Genomics**")
    st.dataframe(gdf_filtered.head(top_n))
    st.write("**Transcriptomics**")
    st.dataframe(tdf_filtered.head(top_n))
    st.write("**Proteomics**")
    st.dataframe(pdf_filtered.head(top_n))

    # Integration
    union_genes = set(gdf_filtered['Gene']) | set(tdf_filtered['Gene'])
    unique_uniprot_ids = extract_uniprot_ids(pdf_filtered['Protein'])
    uniprot_gene_map = map_uniprot_to_gene(unique_uniprot_ids)

    expanded_rows = []
    for _, row in pdf_filtered.iterrows():
        for pid in str(row['Protein']).split(';'):
            pid = pid.strip()
            gene = uniprot_gene_map.get(pid)
            if gene:
                expanded_rows.append({'Protein': pid, 'GeneName': gene})
    expanded_protein_df = pd.DataFrame(expanded_rows)
    protein_gene_map = dict(zip(expanded_protein_df['Protein'], expanded_protein_df['GeneName']))

    # Merge and cluster
    merged_df = pd.merge(gdf_filtered, tdf_filtered, on='Gene')
    merged_df = pd.merge(merged_df, pdf_filtered, on='Gene')
    scaler = StandardScaler()
    features = scaler.fit_transform(merged_df[['CADD', 'logFC', 'AveExpr', 'B', 'Intensity']])
    reducer = umap.UMAP(n_components=2, random_state=42)
    umap_coords = reducer.fit_transform(features)

    kmeans = KMeans(n_clusters=n_clusters, random_state=42)
    clusters = kmeans.fit_predict(umap_coords)
    merged_df['Cluster'] = clusters
    merged_df['UMAP1'], merged_df['UMAP2'] = umap_coords[:, 0], umap_coords[:, 1]

    # Enrichment per cluster
    enrichment_dict = {}
    all_enrichments = []
    libraries = {"Pathway": "Reactome_2016", "Metabolite": "HMDB_Metabolites", "Disease": "DisGeNET"}

    gene_cluster_map = {}
    for c in sorted(merged_df['Cluster'].unique()):
        cluster_genes = merged_df[merged_df['Cluster'] == c]['Gene'].dropna().unique().tolist()
        enrichments = {"Cluster": c}
        for label, lib in libraries.items():
            try:
                enr = enrichr(gene_list=cluster_genes, gene_sets=lib, outdir=None)
                if not enr.results.empty:
                    top_terms = "; ".join(enr.results.sort_values('P-value').head(3)['Term'])
                    enrichments[label] = top_terms
                    for term in enr.results['Term'].head(3):
                        all_enrichments.append((label, term))
            except:
                enrichments[label] = ""
        enrichment_dict[c] = enrichments

        for gene in cluster_genes:
            gene_cluster_map.setdefault((c, gene), {"Protein": [], "Pathway": [], "Metabolite": [], "Disease": []})
            for prot, gname in protein_gene_map.items():
                if gname == gene:
                    gene_cluster_map[(c, gene)]["Protein"].append(prot)
            for label in ["Pathway", "Metabolite", "Disease"]:
                for term in enrichments[label].split(';'):
                    term = term.strip()
                    if term:
                        gene_cluster_map[(c, gene)][label].append(term)

      # Enrichment Summary by Cluster (modified to separate columns for Metabolite, Pathway, Disease)
    # -----------------------------
    st.subheader("📚 Enrichment Summary by Cluster")
    enrichment_summary = []
    for c, enrichments in enrichment_dict.items():
        row = {
            'Cluster': c,
            'Metabolite': enrichments['Metabolite'] if enrichments['Metabolite'] else '',
            'Pathway': enrichments['Pathway'] if enrichments['Pathway'] else '',
            'Disease': enrichments['Disease'] if enrichments['Disease'] else ''
        }
        enrichment_summary.append(row)
    enrichment_summary_df = pd.DataFrame(enrichment_summary)
    st.dataframe(enrichment_summary_df, use_container_width=True)

    # Final barplot for repeated terms
    st.subheader("📊 Frequent Enrichment Terms Across Clusters")
    enrichment_df_all = pd.DataFrame(all_enrichments, columns=['Type', 'Term'])
    bar_data = enrichment_df_all.groupby(['Type', 'Term']).size().reset_index(name='Count')
    fig_bar = px.bar(bar_data, x='Term', y='Count', color='Type', title="Enrichment Term Frequency")
    st.plotly_chart(fig_bar)

    # UMAP with circular cluster boundaries
    st.subheader("🔬 UMAP Visualization with Clusters")
    fig_umap = px.scatter(
        merged_df,
        x='UMAP1', y='UMAP2',
        color=merged_df['Cluster'].astype(str),
        hover_name='Gene',
        hover_data=['CADD', 'logFC', 'Intensity'],
        color_discrete_sequence=px.colors.qualitative.Bold
    )
    for cluster_id in merged_df['Cluster'].unique():
        cluster_data = merged_df[merged_df['Cluster'] == cluster_id]
        x_center, y_center = cluster_data[['UMAP1', 'UMAP2']].mean()
        radius = np.linalg.norm(cluster_data[['UMAP1', 'UMAP2']] - [x_center, y_center], axis=1).max()
        fig_umap.add_shape(type="circle", xref="x", yref="y",
                           x0=x_center - radius, y0=y_center - radius,
                           x1=x_center + radius, y1=y_center + radius,
                           line=dict(color="black", width=2, dash="dash"),
                           fillcolor="rgba(173,216,230,0.1)")
        fig_umap.add_annotation(x=x_center, y=y_center, text=f"Cluster {cluster_id}", showarrow=False, font=dict(size=14))
    st.plotly_chart(fig_umap, use_container_width=True)

    # Network
    st.subheader("🧠 Cluster-Based Interaction Network")
    net = Network(height='800px', width='100%', directed=False)
    net.force_atlas_2based()

    color_map = {"Gene": "gray", "Protein": "gold", "Pathway": "skyblue", "Metabolite": "lightgreen", "Disease": "lightcoral"}
    for (c, gene), info in gene_cluster_map.items():
        net.add_node(gene, label=gene, color=color_map['Gene'])
        for prot in info['Protein']:
            net.add_node(prot, label=prot, color=color_map['Protein'])
            net.add_edge(gene, prot, color='gray')
        for etype in ['Pathway', 'Metabolite', 'Disease']:
            for term in info[etype]:
                net.add_node(term, label=term, color=color_map[etype])
                net.add_edge(gene, term, color='gray')

    with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp_file:
        net.save_graph(tmp_file.name)
        st.components.v1.html(open(tmp_file.name, 'r', encoding='utf-8').read(), height=800)

    st.subheader("🗂️ Network Legend")
    for k, v in color_map.items():
        st.markdown(f"<span style='display:inline-block; width:20px; height:20px; background-color:{v}; margin-right:10px; border-radius:50%;'></span> **{k}**", unsafe_allow_html=True)

    # Final Interaction Table
    st.subheader("📄 Interaction Table from Network")
    rows = []
    for (c, gene), info in gene_cluster_map.items():
        row = {
            'Cluster': c,
            'Gene': gene,
            'Protein': "; ".join(info['Protein']),
            'Metabolite': "; ".join(info['Metabolite']),
            'Pathway': "; ".join(info['Pathway']),
            'Disease': "; ".join(info['Disease'])
        }
        rows.append(row)
    st.dataframe(pd.DataFrame(rows))
