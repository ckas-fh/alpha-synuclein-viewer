import streamlit as st
import py3Dmol
import requests
import streamlit.components.v1 as components

# Page config
st.set_page_config(
    page_title="α-Synuclein 3D Viewer",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 α-Synuclein 3D Protein Viewer")
st.markdown("Interactive visualization of α-synuclein structures and Parkinson's disease mutations")

# Sidebar for controls
st.sidebar.header("Structure Controls")

# PDB structure options for α-synuclein
pdb_options = {
    "1XQ8 - α-Synuclein Fibril": "1XQ8",
    "6H6B - α-Synuclein Monomer": "6H6B", 
    "6CU7 - α-Synuclein in Membrane": "6CU7"
}

selected_structure = st.sidebar.selectbox(
    "Choose α-synuclein structure:",
    options=list(pdb_options.keys())
)
pdb_id = pdb_options[selected_structure]

# Visualization style options
style_options = {
    "Cartoon": "cartoon",
    "Stick": "stick", 
    "Sphere": "sphere",
    "Surface": "surface"
}

selected_style = st.sidebar.selectbox(
    "Visualization style:",
    options=list(style_options.keys()),
    index=0
)
style = style_options[selected_style]

# Color scheme options
color_options = {
    "Spectrum (Rainbow)": "spectrum",
    "Secondary Structure": "sstruc",
    "Chain": "chain",
    "Residue Type": "resn"
}

selected_color = st.sidebar.selectbox(
    "Color scheme:",
    options=list(color_options.keys()),
    index=0
)
color_scheme = color_options[selected_color]

# Known Parkinson's mutations
st.sidebar.markdown("### 🔴 Parkinson's Mutations")
show_mutations = st.sidebar.checkbox("Highlight disease mutations", value=True)

mutations_info = {
    "A53T": {"position": 53, "description": "Most common familial mutation"},
    "A30P": {"position": 30, "description": "Associated with early onset"},
    "E46K": {"position": 46, "description": "Linked to dementia"},
    "H50Q": {"position": 50, "description": "Rare familial variant"},
    "G51D": {"position": 51, "description": "Recently discovered"}
}

if show_mutations:
    selected_mutations = st.sidebar.multiselect(
        "Select mutations to highlight:",
        options=list(mutations_info.keys()),
        default=["A53T", "A30P", "E46K"]
    )

# Create two columns for layout
col1, col2 = st.columns([2, 1])

with col1:
    st.subheader(f"Structure: {selected_structure}")
    
    # Function to create 3D viewer
    def create_3d_viewer(pdb_id, style, color_scheme, mutations=None):
        try:
            # Fetch PDB data
            pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
            response = requests.get(pdb_url)
            pdb_data = response.text
            
            # Create HTML with 3Dmol.js
            mutations_list = mutations if mutations and show_mutations else []
            
            viewer_html = f"""
            <div id="container-01" style="height: 600px; width: 100%; position: relative;"></div>
            
            <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.4/3Dmol-min.js"></script>
            <script>
            let element = document.getElementById('container-01');
            let config = {{ backgroundColor: 'black' }};
            let viewer = $3Dmol.createViewer(element, config);
            
            // PDB data
            let pdbData = `{pdb_data}`;
            
            // Add model
            viewer.addModel(pdbData, "pdb");
            
            // Set style based on selection
            let styleType = "{style}";
            let colorScheme = "{color_scheme}";
            
            if (styleType === "cartoon") {{
                viewer.setStyle({{}}, {{cartoon: {{colorscheme: colorScheme}}}});
            }} else if (styleType === "stick") {{
                viewer.setStyle({{}}, {{stick: {{colorscheme: colorScheme}}}});
            }} else if (styleType === "sphere") {{
                viewer.setStyle({{}}, {{sphere: {{colorscheme: colorScheme}}}});
            }} else if (styleType === "surface") {{
                viewer.setStyle({{}}, {{surface: {{colorscheme: colorScheme, opacity: 0.8}}}});
            }}
            
            // Highlight mutations
            let mutations = {mutations_list};
            let mutationInfo = {mutations_info};
            
            mutations.forEach(function(mutation) {{
                let pos = mutationInfo[mutation].position;
                viewer.addStyle({{resi: pos}}, {{sphere: {{color: 'red', radius: 2.0}}}});
                viewer.addLabel(mutation, {{position: {{resi: pos}}, 
                                         backgroundColor: 'red', 
                                         fontColor: 'white',
                                         fontSize: 12,
                                         borderThickness: 1}});
            }});
            
            // Final setup
            viewer.zoomTo();
            viewer.spin(true);
            viewer.render();
            </script>
            """
            
            return viewer_html
            
        except Exception as e:
            return f"""
            <div style="height: 600px; width: 100%; display: flex; align-items: center; justify-content: center; border: 2px dashed #ccc;">
                <div style="text-align: center;">
                    <h3>Error loading structure</h3>
                    <p>Could not fetch PDB data for {pdb_id}</p>
                    <p>Error: {str(e)}</p>
                    <p>Check your internet connection or try a different structure.</p>
                </div>
            </div>
            """
    
    # Generate and display the viewer
    try:
        mutations_to_highlight = selected_mutations if show_mutations and 'selected_mutations' in locals() else []
        viewer_html = create_3d_viewer(pdb_id, style, color_scheme, mutations_to_highlight)
        components.html(viewer_html, height=650)
        
    except Exception as e:
        st.error(f"Error creating viewer: {e}")
        st.info("There was an issue creating the 3D viewer. Please check your internet connection.")

with col2:
    st.subheader("Protein Information")
    
    # Structure info
    st.markdown(f"**PDB ID:** {pdb_id}")
    st.markdown(f"**Protein:** α-Synuclein")
    st.markdown(f"**Organism:** Homo sapiens")
    st.markdown(f"**Length:** ~140 amino acids")
    
    # Mutation information
    if show_mutations and 'selected_mutations' in locals():
        st.subheader("🔴 Selected Mutations")
        for mutation in selected_mutations:
            with st.expander(f"{mutation} (Position {mutations_info[mutation]['position']})"):
                st.write(mutations_info[mutation]["description"])
                st.write(f"**Position:** {mutations_info[mutation]['position']}")
                st.write("**Impact:** Associated with familial Parkinson's disease")
    
    # Additional info
    st.subheader("About α-Synuclein")
    st.markdown("""
    α-Synuclein is a 140-amino acid protein that:
    - Aggregates to form Lewy bodies in Parkinson's disease
    - Is normally involved in synaptic vesicle trafficking
    - Misfolds and forms toxic oligomers and fibrils
    - Contains three main regions: N-terminal, NAC, C-terminal
    """)
    
    st.subheader("Disease Relevance")
    st.markdown("""
    **Parkinson's Connection:**
    - Mutations cause familial Parkinson's disease
    - Wild-type protein also aggregates in sporadic cases
    - Target for therapeutic interventions
    - Key protein in Lewy body formation
    """)

# Footer with next steps
st.markdown("---")
st.subheader("🚀 Next Development Steps")
col1, col2, col3 = st.columns(3)

with col1:
    st.markdown("""
    **Phase 1: Core Viewer**
    - ✅ Structure selection
    - ✅ Basic styling options
    - ✅ Mutation highlighting
    - 🔄 Interactive 3D viewer
    """)

with col2:
    st.markdown("""
    **Phase 2: Advanced Features**
    - 📋 Binding site mapping
    - 📊 Property calculations
    - 🔄 Structure comparison
    - 💾 Save/export options
    """)

with col3:
    st.markdown("""
    **Phase 3: Research Tools**
    - 📚 Literature integration
    - 🧪 Drug interaction sites
    - 📈 Aggregation predictions
    - 👥 User accounts & sharing
    """)

st.info("💡 **To run this locally:** `pip install streamlit py3Dmol requests` then `streamlit run app.py`")