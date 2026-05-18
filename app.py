import streamlit as st
import pandas as pd
from pyvis.network import Network
import streamlit.components.v1 as components
import networkx as nx
import requests
import time
import os

st.set_page_config(page_title="Stress-Path: Final Build", layout="wide")

st.title("🌱 Stress-Path: Metabolic Mapper")
st.markdown("### Upload your Stress Dataset (CSV)")


@st.cache_data(show_spinner=False)
def get_kegg_id(compound_name):
    query = str(compound_name).strip().lower()
    url = f"http://rest.kegg.jp/find/compound/{query}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        if response.text:
            return response.text.split('\n')[0].split('\t')[0].replace('cpd:', '')
    except requests.exceptions.RequestException:
        return None
    return None

@st.cache_data(show_spinner=False)
def get_chickpea_pathways(kegg_id):
    if not kegg_id: return []
    url = f"http://rest.kegg.jp/link/pathway/{kegg_id}"
    pathways = []
    try:
        response = requests.get(url)
        response.raise_for_status()
        if response.text:
            for line in response.text.strip().split('\n'):
                parts = line.split('\t')
                if len(parts) >= 2:
                    pathway_id = parts[1].replace('path:', '')
                    if pathway_id.startswith('car') or pathway_id.startswith('map'):  
                        pathways.append(pathway_id)
    except requests.exceptions.RequestException:
        pass
    time.sleep(0.1) 
    return pathways
# --------------------------------------------------------
# ---------------------------

st.sidebar.header("Settings")
vip_threshold = st.sidebar.slider("VIP Threshold", 0.0, 2.0, 1.0, step=0.1)

uploaded_file = st.file_uploader("Choose a file", type="csv")

if uploaded_file is not None:
    data = pd.read_csv(uploaded_file)
    st.write("### Data Preview")
    st.dataframe(data.head())
    
    significant_genes = data[data['VIP'] > vip_threshold]
    st.success(f"Found {len(significant_genes)} significant metabolic changes!")
    
    csv_buffer = significant_genes.to_csv(index=False).encode('utf-8')
    st.download_button("📥 Download Filtered Stress Responders", data=csv_buffer, file_name="significant_metabolites_report.csv", mime="text/csv")

    st.write("### Significant Metabolic Switches")
    st.dataframe(significant_genes[['Compounds', 'Class I', 'VIP']])

    st.write("### 🧬 Live KEGG Database Integration & Comparative Analysis")
    
    compound_data = {}
    nodes = significant_genes['Compounds'].tolist()
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    with st.spinner("Connecting to KEGG REST API (Kyoto Encyclopedia of Genes and Genomes)..."):
        for i, compound in enumerate(nodes):
            status_text.text(f"Querying KEGG for: {compound}...")
            k_id = get_kegg_id(compound)
            paths = get_chickpea_pathways(k_id)
            
            bio_class = significant_genes[significant_genes['Compounds'] == compound]['Class I'].values[0]
            vip_score = significant_genes[significant_genes['Compounds'] == compound]['VIP'].values[0]
            reg_type = significant_genes[significant_genes['Compounds'] == compound]['Type'].values[0]
            
            compound_data[compound] = {
                "kegg_id": k_id, 
                "pathways": paths, 
                "class": bio_class, 
                "vip": vip_score,
                "type": reg_type
            }
            
            progress_bar.progress((i + 1) / len(nodes))
            
    status_text.text("KEGG Mapping Complete! Rendering Biological Network...")
    time.sleep(1) 
    status_text.empty()
    progress_bar.empty()

    nx_graph = nx.Graph()

    for compound, info in compound_data.items():
        nx_graph.add_node(compound, group=info['class'], vip=info['vip'], kegg_id=info['kegg_id'], type=info['type'])

    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            node_a, node_b = nodes[i], nodes[j]
            paths_a, paths_b = set(compound_data[node_a]["pathways"]), set(compound_data[node_b]["pathways"])
            
            shared_paths = paths_a.intersection(paths_b)
            if shared_paths:
                path_label = list(shared_paths)[0]
                nx_graph.add_edge(node_a, node_b, shared_pathway=path_label)

    pos = nx.spring_layout(nx_graph, seed=42, k=0.8)

    net = Network(height='600px', width='100%', bgcolor='#ffffff', font_color='black')
    net.toggle_physics(False)

    for node in nx_graph.nodes():
        vip_score = nx_graph.nodes[node]['vip']
        bio_class = nx_graph.nodes[node]['group']
        kegg_id = nx_graph.nodes[node]['kegg_id']
        reg_type = str(nx_graph.nodes[node]['type']).strip().lower() # Get 'up' or 'down'
        
        if reg_type == 'up':
            node_color = '#e74c3c'  # Strong Red
        elif reg_type == 'down':
            node_color = '#3498db'  # Strong Blue
        else:
            node_color = '#95a5a6'  # Gray fallback
            
        node_size = int(vip_score * 15)
        
        hover_title = f"Class: {bio_class}\nRegulation: {reg_type.upper()}\nKEGG ID: {kegg_id if kegg_id else 'Unmapped'}"
        
        x_coord = float(pos[node][0] * 800)
        y_coord = float(pos[node][1] * 800)
        
        net.add_node(node, label=node, title=hover_title, color=node_color, size=node_size, x=x_coord, y=y_coord)

    for edge in nx_graph.edges():
        shared_pathway = nx_graph.edges[edge]['shared_pathway']
        net.add_edge(edge[0], edge[1], title=f"Shared KEGG Pathway: {shared_pathway}", color="#a9cce3", width=2)

    path = "html_files"
    if not os.path.exists(path): os.makedirs(path)
    net.save_graph(f"{path}/network.html")
    with open(f"{path}/network.html", 'r', encoding='utf-8') as f:
        components.html(f.read(), height=650)

    st.write("### 🔑 Network Topology Insights")
    st.caption("Calculated via NetworkX Degree Centrality based on KEGG mapping")
    
    if len(nx_graph.edges()) > 0:
        centrality = nx.degree_centrality(nx_graph)
        metrics_df = pd.DataFrame(list(centrality.items()), columns=['Metabolite', 'Centrality Score'])
        st.dataframe(metrics_df.sort_values(by='Centrality Score', ascending=False).head(5))
    else:
        st.warning("No shared KEGG pathways found. This highlights a gap in global database annotation for these specific secondary metabolites.")

else:
    st.info("Awaiting CSV upload. Please upload the metabolite_data.csv file.")