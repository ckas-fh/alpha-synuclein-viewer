import streamlit as st
import requests
import streamlit.components.v1 as components
import numpy as np
from typing import Dict, List, Tuple
import re

# Your AggregationPredictor class
class AggregationPredictor:
    """
    AI-powered protein aggregation risk predictor for α-synuclein.
    Uses sequence-based features to predict aggregation-prone regions.
    """
    
    def __init__(self):
        # Amino acid properties for feature calculation
        self.hydrophobicity = {
            'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
            'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
            'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
        }
        
        self.charge = {
            'A': 0, 'R': 1, 'N': 0, 'D': -1, 'C': 0, 'Q': 0, 'E': -1,
            'G': 0, 'H': 0.5, 'I': 0, 'L': 0, 'K': 1, 'M': 0, 'F': 0,
            'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0
        }
        
        # Beta-sheet propensity (aggregation-prone secondary structure)
        self.beta_propensity = {
            'A': 0.83, 'R': 0.93, 'N': 0.89, 'D': 0.54, 'C': 1.19,
            'Q': 1.10, 'E': 0.37, 'G': 0.75, 'H': 0.87, 'I': 1.60,
            'L': 1.30, 'K': 0.74, 'M': 1.05, 'F': 1.38, 'P': 0.55,
            'S': 0.75, 'T': 1.19, 'W': 1.37, 'Y': 1.47, 'V': 1.70
        }
        
        # Known aggregation-prone motifs in α-synuclein
        self.aggregation_motifs = [
            'NAC',      # Non-amyloid component region
            'KTKE',     # Known aggregation motif
            'GVLY',     # Hydrophobic cluster
            'VLYVG',    # Core aggregation region
            'GAVVT',    # Another aggregation-prone sequence
        ]
        
        # α-synuclein sequence (human, UniProt P37840)
        self.alpha_syn_sequence = (
            "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQV"
            "TNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYE"
            "MPSEEGYQDYEPEA"
        )
    
    def extract_features(self, sequence: str) -> Dict[str, List[float]]:
        """Extract aggregation-relevant features from protein sequence."""
        features = {
            'hydrophobicity': [],
            'charge': [],
            'beta_propensity': [],
            'local_hydrophobicity': [],
            'charge_clustering': [],
            'motif_score': []
        }
        
        sequence = sequence.upper()
        
        for i, aa in enumerate(sequence):
            # Basic properties
            features['hydrophobicity'].append(self.hydrophobicity.get(aa, 0))
            features['charge'].append(self.charge.get(aa, 0))
            features['beta_propensity'].append(self.beta_propensity.get(aa, 0))
            
            # Local hydrophobicity (sliding window average)
            window = 5
            start = max(0, i - window//2)
            end = min(len(sequence), i + window//2 + 1)
            local_hydro = np.mean([self.hydrophobicity.get(sequence[j], 0) 
                                 for j in range(start, end)])
            features['local_hydrophobicity'].append(local_hydro)
            
            # Charge clustering (local charge concentration)
            local_charge = np.sum([abs(self.charge.get(sequence[j], 0)) 
                                 for j in range(start, end)])
            features['charge_clustering'].append(local_charge)
            
            # Motif score (proximity to known aggregation motifs)
            motif_score = 0
            for motif in self.aggregation_motifs:
                for j in range(len(sequence) - len(motif) + 1):
                    if sequence[j:j+len(motif)] == motif:
                        distance = abs(i - (j + len(motif)//2))
                        motif_score += max(0, 5 - distance) / 5  # Decay with distance
            features['motif_score'].append(motif_score)
        
        return features
    
    def calculate_aggregation_risk(self, sequence: str) -> List[float]:
        """Calculate aggregation risk score for each position (0-100)."""
        features = self.extract_features(sequence)
        risk_scores = []
        
        for i in range(len(sequence)):
            # Weighted combination of features
            risk = 0
            
            # High hydrophobicity increases risk
            hydro_risk = max(0, features['hydrophobicity'][i]) * 10
            
            # High beta-sheet propensity increases risk
            beta_risk = features['beta_propensity'][i] * 15
            
            # Local hydrophobic clusters increase risk
            cluster_risk = max(0, features['local_hydrophobicity'][i]) * 8
            
            # Low charge clustering increases risk (neutral regions aggregate more)
            charge_risk = max(0, 3 - features['charge_clustering'][i]) * 5
            
            # Proximity to known motifs increases risk
            motif_risk = features['motif_score'][i] * 25
            
            # Combine all risk factors
            total_risk = hydro_risk + beta_risk + cluster_risk + charge_risk + motif_risk
            
            # Normalize to 0-100 scale
            risk_scores.append(min(100, max(0, total_risk)))
        
        return risk_scores
    
    def get_high_risk_regions(self, sequence: str, threshold: float = 60) -> List[Tuple[int, int]]:
        """Identify continuous high-risk regions above threshold."""
        risk_scores = self.calculate_aggregation_risk(sequence)
        regions = []
        start = None
        
        for i, score in enumerate(risk_scores):
            if score >= threshold and start is None:
                start = i
            elif score < threshold and start is not None:
                regions.append((start, i-1))
                start = None
        
        # Handle case where high-risk region extends to end
        if start is not None:
            regions.append((start, len(sequence)-1))
        
        return regions
    
    def get_color_map(self, risk_scores: List[float]) -> List[str]:
        """Convert risk scores to color map for visualization."""
        colors = []
        for score in risk_scores:
            if score < 20:
                colors.append('#00FF00')  # Green - low risk
            elif score < 40:
                colors.append('#80FF00')  # Yellow-green
            elif score < 60:
                colors.append('#FFFF00')  # Yellow
            elif score < 80:
                colors.append('#FF8000')  # Orange
            else:
                colors.append('#FF0000')  # Red - high risk
        return colors

# Page config
st.set_page_config(
    page_title="α-Synuclein 3D Viewer",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 α-Synuclein 3D Protein Viewer")
st.markdown("Interactive visualization of α-synuclein structures and Parkinson's disease mutations with AI aggregation prediction")

# Initialize aggregation predictor
@st.cache_resource
def load_predictor():
    return AggregationPredictor()

predictor = load_predictor()

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
    "Aggregation Risk (AI)": "aggregation",
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

# AI Aggregation Analysis
st.sidebar.markdown("### 🤖 AI Aggregation Analysis")
show_aggregation = st.sidebar.checkbox("Show aggregation risk prediction", value=True)

if show_aggregation:
    risk_threshold = st.sidebar.slider(
        "Risk threshold for highlighting:",
        min_value=30, max_value=90, value=60, step=10
    )
    
    # Calculate aggregation risk
    sequence = predictor.alpha_syn_sequence
    risk_scores = predictor.calculate_aggregation_risk(sequence)
    high_risk_regions = predictor.get_high_risk_regions(sequence, risk_threshold)
    avg_risk = np.mean(risk_scores)
    max_risk = max(risk_scores)

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
            
            # Generate aggregation coloring if selected
            aggregation_coloring = ""
            if color_scheme == "aggregation" and show_aggregation:
                risk_scores = predictor.calculate_aggregation_risk(predictor.alpha_syn_sequence)
                colors = predictor.get_color_map(risk_scores)
                for i, color in enumerate(colors):
                    aggregation_coloring += f"viewer.addStyle({{resi: {i+1}}}, {{cartoon: {{color: '{color}'}}}});\n"
            
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
            
            if (colorScheme === "aggregation") {{
                // Apply AI aggregation coloring
                viewer.setStyle({{}}, {{{style}: {{}}}});
                {aggregation_coloring}
            }} else {{
                // Apply standard coloring
                if (styleType === "cartoon") {{
                    viewer.setStyle({{}}, {{cartoon: {{colorscheme: colorScheme}}}});
                }} else if (styleType === "stick") {{
                    viewer.setStyle({{}}, {{stick: {{colorscheme: colorScheme}}}});
                }} else if (styleType === "sphere") {{
                    viewer.setStyle({{}}, {{sphere: {{colorscheme: colorScheme}}}});
                }} else if (styleType === "surface") {{
                    viewer.setStyle({{}}, {{surface: {{colorscheme: colorScheme, opacity: 0.8}}}});
                }}
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
    
    # AI Aggregation Analysis Results
    if show_aggregation:
        st.subheader("🤖 AI Aggregation Analysis")
        
        col_a, col_b = st.columns(2)
        with col_a:
            st.metric("Average Risk", f"{avg_risk:.1f}/100")
        with col_b:
            st.metric("Max Risk", f"{max_risk:.1f}/100")
        
        # High-risk regions
        if high_risk_regions:
            st.markdown("**🔥 High-Risk Regions:**")
            for start, end in high_risk_regions:
                region_seq = sequence[start:end+1]
                avg_region_risk = np.mean(risk_scores[start:end+1])
                st.markdown(f"• **Positions {start+1}-{end+1}:** {region_seq}")
                st.markdown(f"  Risk: {avg_region_risk:.1f}/100")
        else:
            st.markdown("✅ **No high-risk regions found** above threshold")
        
        # Risk interpretation
        with st.expander("🧠 Understanding Aggregation Risk"):
            st.markdown("""
            **Risk Score Interpretation:**
            - **0-20:** Low risk (green) - stable regions
            - **20-40:** Moderate risk (yellow-green)
            - **40-60:** Elevated risk (yellow)  
            - **60-80:** High risk (orange) - prone to aggregation
            - **80-100:** Very high risk (red) - likely aggregation sites
            
            **AI Features Used:**
            - Amino acid hydrophobicity
            - Beta-sheet formation propensity
            - Local hydrophobic clustering
            - Charge distribution analysis
            - Known aggregation motif proximity
            """)
    
    # Mutation information
    if show_mutations and 'selected_mutations' in locals():
        st.subheader("🔴 Selected Mutations")
        for mutation in selected_mutations:
            with st.expander(f"{mutation} (Position {mutations_info[mutation]['position']})"):
                st.write(mutations_info[mutation]["description"])
                st.write(f"**Position:** {mutations_info[mutation]['position']}")
                st.write("**Impact:** Associated with familial Parkinson's disease")
                
                # Show aggregation risk for this position if available
                if show_aggregation:
                    pos = mutations_info[mutation]['position'] - 1  # Convert to 0-indexed
                    if pos < len(risk_scores):
                        risk_at_pos = risk_scores[pos]
                        st.write(f"**AI Aggregation Risk:** {risk_at_pos:.1f}/100")
                        if risk_at_pos > 60:
                            st.warning("⚠️ This mutation is in a high-risk aggregation region!")
                        else:
                            st.info("ℹ️ This mutation is in a lower-risk region.")
    
    # Additional info
    st.subheader("About α-Synuclein")
    st.markdown("""
    α-Synuclein is a 140-amino acid protein that:
    - Aggregates to form Lewy bodies in Parkinson's disease
    - Is normally involved in synaptic vesicle trafficking
    - Misfolds and forms toxic oligomers and fibrils
    - Contains three main regions: N-terminal, NAC, C-terminal
    """)

# Footer
st.markdown("---")
st.subheader("🚀 AI-Powered Features")
col1, col2, col3 = st.columns(3)

with col1:
    st.markdown("""
    **✅ Core AI Features**
    - 🧠 Advanced aggregation prediction
    - 🎨 Risk-based protein coloring
    - 📊 Multi-feature analysis
    - 🔍 High-risk region detection
    """)

with col2:
    st.markdown("""
    **🔬 Scientific Features**
    - 🧬 3D structure visualization
    - 🔴 Disease mutation mapping
    - 📈 Quantitative risk scoring
    - 💡 Interpretable predictions
    """)

with col3:
    st.markdown("""
    **💻 Technical Features**
    - ⚡ Real-time analysis
    - 🎛️ Interactive controls
    - 📱 Responsive design
    - 🌐 Web-based deployment
    """)

st.info("💡 **Portfolio Highlight:** This project demonstrates AI/ML application to real biological problems - exactly what biotech companies are looking for!")