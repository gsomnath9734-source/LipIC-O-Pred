from sklearn.linear_model import LinearRegression
import streamlit as st
from streamlit_ketcher import st_ketcher
import numpy as np
import pandas as pd
import base64
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from rasar import ra_pred
from lib.leverage import leverage_calculator
import os

qsar_data = pd.read_excel("lib/qsar.xlsx", index_col=0)
rasar_data = pd.read_excel("lib/rasar.xlsx", index_col=0)

mlr_model = LinearRegression().fit(rasar_data.iloc[:, :-1], rasar_data.iloc[:, -1])

selected_descriptors = [
    "EState_VSA3",
    "FractionCSP3",
    "BalabanJ",
    "SMR_VSA1",
    "NumSaturatedHeterocycles",
    "SlogP_VSA11",
    "EState_VSA5",
    "PEOE_VSA3"
]

desc_dict = dict(Descriptors.descList)
selected_funcs = [(name, desc_dict[name]) for name in selected_descriptors]

def smiles_to_descriptor_df(smiles_list):
    data = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            row = {name: None for name, _ in selected_funcs}
        else:
            row = {name: func(mol) for name, func in selected_funcs}
        data.append(row)
    return pd.DataFrame(data, index=smiles_list)

st.set_page_config(
    page_title="Lipase Inhibitory Concentration for Obesity Predictor",
    layout="wide"
)

import requests

def get_base64_image_from_url(url):
    response = requests.get(url)
    response.raise_for_status()
    return base64.b64encode(response.content).decode()

img_base64 = get_base64_image_from_url(
    "https://raw.githubusercontent.com/gsomnath9734-source/LipIC-O-Pred/main/Logo.png"
)

st.markdown(
    f"""
    <div style="text-align: center;">
        <img src="data:image/png;base64,{img_base64}" style="width: 400px; height: auto;">
    </div>
    """,
    unsafe_allow_html=True
)

st.markdown(
    "<h4 style='text-align: center;'> Lipase Inhibitory Concentration for Obesity Predictor </h4>",
    unsafe_allow_html=True
)

with st.expander(
    "**LipIC-O-Pred v1.0** is a user-friendly, regression-based pIC₅₀ prediction tool for small molecules against pancreatic lipase in obesity management, using RA function (LK) RASAR descriptor that integrates chemical similarity information of the close source compounds using the Laplacian Kernel approach.",
    expanded=True
):
    st.markdown(
        "<p style='text-align: center;'>"
        "Name, class and definitions of QSAR descriptors used to calculate RA function (LK) in LipIC-O-Pred v1.0"
        "</p>",
        unsafe_allow_html=True
    )

    st.markdown("""
| **Descriptor** | **Class** | **Definition** |
|---------------|-----------|----------------|
| EState_VSA3 | Estate VSA | The sum of Van der Waals surface area contributions for atoms whose electrotopological state values (x) fall within the range of 0.29 ≤ x < 0.72. |
| FractionCSP3 | Molecular property | Fraction of sp³ hybridized carbon atoms. |
| BalabanJ | Topographical | The sum of the reciprocal square roots of the distance sums for each edge, normalized by the cyclomatic number. |
| SMR_VSA1 | MOE-type | The sum of Van der Waals surface area contributions for atoms with molar refractivity (x) in the range of −∞ < x < 1.29. |
| NumSaturatedHeterocycles | Topographical | Number of saturated heterocycles in a molecule. |
| SlogP_VSA11 | VSA descriptor | The sum of Van der Waals surface area contributions of atoms with Wildman–Crippen lipophilicity (x) in the range of 0.50 ≤ x < 0.60. |
| EState_VSA5 | Estate VSA | The sum of Van der Waals surface area for atoms with electrotopological state values (x) in the range of 1.17 ≤ x < 1.54. |
| PEOE_VSA3 | MOE-type | Partial Equalized Orbital Electronegativity VSA with partial charge (x) in the range of −0.25 ≤ x < −0.20. |
""")
    
st.space()

def load_file(uploaded_data):
    if uploaded_data is not None:
        return pd.read_excel(uploaded_data, index_col=0)
    return None

def main_function():
    if smiles_input and Chem.MolFromSmiles(smiles_input) is not None and target_df is not None:
        with col2:
            st.warning("Please provide input by drawing a query molecule or uploading a batch Excel file to proceed with the analysis.")
    elif smiles_input and Chem.MolFromSmiles(smiles_input) is not None:
        desc_df = smiles_to_descriptor_df([smiles_input])
        rasar_desc = ra_pred(df1=qsar_data, df2 = desc_df).weighted_prediction(method="Laplacian Kernel", ctc=4, gamma=0.25)
        rasar_desc.name = "RA function(LK)"
        predicted_pic50 = mlr_model.predict(rasar_desc.values.reshape(1, -1))
       
        __, ad_status, __ = leverage_calculator(data1=rasar_data.iloc[:, :-1], data2=rasar_desc.to_frame().T)
        
        with col2:
            st.info("Analysis completed successfully! Find the results below.")
            mol = Chem.MolFromSmiles(smiles_input)
           
            img = Draw.MolToImage(mol, size=(300, 300))
            left, center, right = st.columns([1, 2, 1])
            with center:
                st.image(img, width=300)
            st.markdown(f"**Predicted pIC₅₀ Value:** {predicted_pic50[0]:.2f}")
            if predicted_pic50[0] >=4.7599:
                st.markdown("##The compound is active based on the training set mean")
            else:
                st.markdown("##The compound is inactive based on the training set mean")
            st.markdown(f"**Applicability Domain (AD) Status:** {ad_status['AD Status'].values[0]}")
            st.markdown(
            """
            <div style="font-size:16px; text-align: justify;">
            <b>Note:</b> The predicted pIC₅₀ (μM) and AD status (leverage method) are obtained using the q-RASAR model (RA function descriptor). 
            “Within AD” indicates reliable predictions within the model’s chemical space, while “Outside AD” suggests the result should be interpreted with caution due to limited representation in the training data.
            </div>
            """,
            unsafe_allow_html=True
            )
            
    elif target_df is not None:
        desc_df = smiles_to_descriptor_df(target_df.iloc[:, 0].tolist())
        rasar_desc = ra_pred(df1=qsar_data, df2 = desc_df).weighted_prediction(method="Laplacian Kernel", ctc=4, gamma=0.25)
        rasar_desc.name = "RA function(LK)"
        rasar_desc_df = rasar_desc.to_frame()
        print(rasar_desc_df)
        print(desc_df)
        predicted_pic50 = mlr_model.predict(rasar_desc.values.reshape(len(rasar_desc), -1))
        __, ad_status, __ = leverage_calculator(data1=rasar_data.iloc[:, :-1], data2=rasar_desc_df[["RA function(LK)"]])
        result_df = pd.DataFrame({
            "Sl. No.": target_df.index.values,
            "Query SMILES": target_df.iloc[:, 0],
            "Predicted pIC₅₀ Value": predicted_pic50,
            "Predicted Class": np.where(predicted_pic50 >= 4.7599, "active", "inactive"),
            "AD Status": ad_status["AD Status"].values
        })
        with col2:
            st.info("Analysis completed successfully! Find the results below.")
        with col2:
            st.dataframe(result_df, hide_index=True)
            st.markdown(
            """
            <div style="font-size:16px; text-align: justify;">
            <b>Note:</b> The predicted pIC₅₀ (μM) and AD status (leverage method) are obtained using the q-RASAR model (RA function descriptor). 
            “Within AD” indicates reliable predictions within the model’s chemical space, while “Outside AD” suggests the result should be interpreted with caution due to limited representation in the training data.
            </div>
            """,
            unsafe_allow_html=True
            )
    else:
        with col2:
            st.warning("Please provide input by drawing a query molecule or uploading a batch Excel file to proceed with the analysis.")


col1, col2 = st.columns(2)

with col1:
    st.markdown("##### Draw Query Molecule and Apply")
    smi_code = st_ketcher()
    input_label = "Enter or Edit SMILES and Run Analysis"
    smiles_input = st.text_input(input_label, value=smi_code)
    

    st.markdown("<h2 style='text-align: center;'>Or</h2>", unsafe_allow_html=True)
    
    
    input_label1 = "Upload SMILES Input (.xlsx) and Run Analysis"
    st.markdown(f"<p style='font-size: 18px; font-weight: 600; color: #31333f;'>{input_label1}</p>", unsafe_allow_html=True)
    target_file = st.file_uploader("", type=["xlsx"])
    target_df = load_file(target_file)

    subcol1, subcol2, subcol3 = st.columns([1, 1, 1])

    with subcol2:
        st.button("🧮 Run Analysis", use_container_width=True, on_click=main_function)
 
with col2:
    
    subcol1, subcol2, subcol3 = st.columns([1, 2, 1])

    with subcol2:
        st.markdown("#### Result(s) of Analysis")


st.markdown("""
<style>
/* ── Google Font ── */
@import url('https://fonts.googleapis.com/css2?family=DM+Serif+Display&family=DM+Sans:wght@400;500;600&display=swap');
 
/* ── Panel wrapper ── */
.panel-box {
    background: linear-gradient(135deg, #f8fafc 0%, #eef2ff 100%);
    border: 1.5px solid #c7d2fe;
    border-radius: 16px;
    padding: 28px 32px;
    margin-top: 12px;
    font-family: 'DM Sans', sans-serif;
    box-shadow: 0 4px 24px rgba(99,102,241,0.08);
    animation: fadeSlide 0.3s ease;
}
.panel-box.contact {
    background: linear-gradient(135deg, #fefce8 0%, #fef9c3 100%);
    border-color: #fcd34d;
    box-shadow: 0 4px 24px rgba(251,191,36,0.10);
}
@keyframes fadeSlide {
    from { opacity: 0; transform: translateY(-10px); }
    to   { opacity: 1; transform: translateY(0); }
}
 
/* ── Panel title ── */
.panel-title {
    font-family: 'DM Serif Display', serif;
    font-size: 1.55rem;
    color: #3730a3;
    margin: 0 0 20px 0;
    border-bottom: 2px solid #c7d2fe;
    padding-bottom: 10px;
    letter-spacing: -0.3px;
}
.panel-title.contact { color: #92400e; border-color: #fcd34d; }
 
/* ── Section heading ── */
.sec-heading {
    font-size: 0.78rem;
    font-weight: 600;
    letter-spacing: 0.12em;
    text-transform: uppercase;
    color: #6366f1;
    margin: 22px 0 8px 0;
}
.sec-heading.contact { color: #b45309; }
 
/* ── Body text ── */
.body-text {
    font-size: 0.93rem;
    color: #374151;
    line-height: 1.7;
    margin: 0;
}
 
/* ── Step cards ── */
.step-row {
    display: flex;
    gap: 14px;
    align-items: flex-start;
    margin-bottom: 10px;
}
.step-label {
    background: #6366f1;
    color: white;
    font-size: 0.78rem;
    font-weight: 700;
    border-radius: 20px;
    padding: 3px 12px;
    white-space: nowrap;
    margin-top: 2px;
    letter-spacing: 0.04em;
}
.step-desc {
    font-size: 0.91rem;
    color: #374151;
    line-height: 1.6;
}
 
/* ── Contact card ── */
.contact-card {
    background: white;
    border: 1px solid #fcd34d;
    border-radius: 10px;
    padding: 14px 18px;
    margin-bottom: 12px;
    font-size: 0.91rem;
    color: #374151;
    line-height: 1.7;
}
.contact-card strong { color: #92400e; }
</style>
""", unsafe_allow_html=True)
   

if "show_manual" not in st.session_state:
    st.session_state.show_manual = False
if "show_info" not in st.session_state:
    st.session_state.show_info = False
 
st.divider()
 
col11, col12, col13 = st.columns([1, 6, 1])
 
with col11:
    if st.button("📖 User Manual", use_container_width=True):
        st.session_state.show_manual = not st.session_state.show_manual
        st.session_state.show_info = False
 
with col13:
    if st.button("ℹ️ Contact Info", use_container_width=True):
        st.session_state.show_info = not st.session_state.show_info
        st.session_state.show_manual = False

if st.session_state.show_manual:
    st.markdown("""
    <div class="panel-box">
 
      <p class="panel-title"> User Manual</p>
 
      <!-- 1. Objective -->
      <p class="sec-heading">Objective</p>
      <p class="body-text">
      The objective of this tool is to develop a user-friendly web-based platform for predicting the pancreatic lipase inhibitory concentration (pIC₅₀ values) of query compounds using a q-RASAR modeling approach. This tool aims to facilitate rapid and reliable screening of bioactive molecules, thereby supporting the identification and optimization of potential anti-obesity agents while reducing the need for extensive experimental studies.
      </p>

      <br>
     
      <!-- 2. System Requirements -->
      <p class="sec-heading">System Requirements</p>
      <p class="body-text">
        <strong>Platform:</strong> Web-based (Runs in any modern browser)<br>
        <strong>Python:</strong> v3.8 or higher<br>
        <strong>Required Packages:</strong> streamlit, streamlit_ketcher, pandas, numpy, rdkit, scikit-learn, pillow, rasar
      </p>

      <br>
        
      <!-- 3. User Guide -->
      <p class="sec-heading">User Guide</p>

      <br>
 
      <div class="step-row">
        <span class="step-label">Step 1</span>
        <span class="step-desc"><strong>Draw or Input Query Molecule</strong> — Begin by drawing the structure of your query compound using the molecular editor provided in the interface. Once the structure is complete, click on the <strong>Apply</strong> button to automatically generate the corresponding SMILES notation. Alternatively, users can directly enter or modify the SMILES string in the input field. After providing the input, click on <strong>Run Analysis</strong> to start the prediction process. The predicted pIC₅₀ value along with other relevant results will be displayed in the panel on the right-hand side of the interface.</span>
      </div>

      <div class="step-row">
        <span class="step-label">Step 2</span>
        <span class="step-desc"><strong>Batch Processing for Multiple Compounds</strong> — For analyzing a large number of compounds simultaneously, users can prepare an Excel file in .xlsx format. The file should contain two columns, where the first column represents the serial number (Sl. No.) and the second column contains the corresponding SMILES strings of the compounds. This file can then be uploaded using the “Upload SMILES Input (.xlsx)” option available in the interface. After uploading, click on <strong>Run Analysis</strong> to process all the compounds at once, and the results will be generated and displayed on the right-side panel.</span>
      </div>

      <div class="step-row">
        <span class="step-label">Step 3</span>
        <span class="step-desc"><strong>View and Download Results</strong> — Once the analysis is complete, the results can be viewed directly within the interface. The output will include important details such as the serial number, query SMILES, predicted pIC₅₀ values, and Applicability Domain (AD) status. Users also have the option to download the results in .csv format for batch processing, which can be used for further analysis, record-keeping, or reporting purposes.</span>
      </div>

      <div class="step-row">
        <span class="step-label">Step 4</span>
        <span class="step-desc"><strong>Interpretation of Results</strong> — The predicted pIC₅₀ value reflects the inhibitory potential of the compound against pancreatic lipase in <strong>micromolar (μM)</strong> unit, where higher values indicate stronger inhibitory activity. The AD status provides information about the reliability of the prediction by indicating whether the query compound falls within the chemical space of the model. Users are encouraged to consider the pIC₅₀ value of <strong>Orlistat</strong> (between 4 to 8 depending on experimental condition variations), an FDA-approved compound and the above parameters while interpreting the results for better decision-making.</span>
      </div>

      <div class="step-row">
        <span class="step-label">Step 5</span>
        <span class="step-desc">
        <strong>Reset and Perform New Analysis</strong> — After completing an analysis, users can click on the <strong>Reset</strong> button to clear the current input or drawn structure. This allows the interface to be reused for a new query. Users can then proceed with drawing a new molecule or uploading another dataset to perform additional predictions.</span>
      </div>

      <br>
 
      <!-- 4. Conditions -->
      <p class="sec-heading">Conditions of Use</p>
      <p class="body-text">
      This software has been developed by the DTC Lab and is intended solely for research purposes. For any inconvenience related to the system or calculations, please feel free to contact the individual listed in the contact information section.</p>
        
    <br>

    <!-- 5. How to Cite -->
   <p class="sec-heading">How to Cite</p>
   <p class="body-text">
   1. A. Banerjee, K. Roy, “First report of q-RASAR modeling toward an approach of easy interpretability and efficient transferability” Mol. Divers. 2022, 26, 2847–2862<br>
   2. S. Pore, K. Roy, “‘intelligent Read Across (iRA)’ - A tool for read-across-based toxicity prediction of nanoparticles” Comput. Struct. Biotechnol. J. 2025, 29, 186–200.
   </p>
    """, unsafe_allow_html=True)

if st.session_state.show_info:
    st.markdown("""
    <div class="panel-box contact">
 
      <p class="panel-title contact"> Contact Information</p>
 
      <p class="sec-heading contact"> Contact </p>
      <div class="contact-card">
        <strong>Prof. (Dr.) Kunal Roy</strong><br>
        <a href="mailto:kunal.roy@jadavpuruniversity.in">kunal.roy@jadavpuruniversity.in</a><br>
        Drug Theoretics and Cheminformatics Laboratory<br>
        Jadavpur University, Kolkata, IN
      </div>
      </p>
                
      <p class="sec-heading contact">Developers</p>
      <div class="contact-card">
        <strong>Somnath Ghosh & Souvik Pore</strong><br>
        <a href="mailto:gsomnath9734@gmail.com">gsomnath9734@gmail.com</a> & 
        <a href="mailto:souvikpore123@gmail.com">souvikpore123@gmail.com</a><br>
        Drug Theoretics and Cheminformatics Laboratory<br>
        Jadavpur University, Kolkata, IN
      </div>
      </p>
    </div>
    """, unsafe_allow_html=True)

    #This platform has been developed by
    #Somnath Ghosh & Souvik Pore
    #DTC laboratory
    #Jadavpur University, Kolkata, India
