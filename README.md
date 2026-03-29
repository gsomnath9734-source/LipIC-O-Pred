# 💊 LipIC-O-Pred

**LipIC-O-Pred** is an interactive, web-based tool for predicting the **pIC50** (negative log of half-maximal inhibitory concentration) of drug-like molecules. Built with [Streamlit](https://streamlit.io/), it provides a seamless cheminformatics workflow — from structure drawing to prediction and applicability domain analysis — all in the browser.

---

## ✨ Features

- 🧪 **Interactive Molecule Input** — Draw molecules directly in the browser using the integrated Ketcher structure editor (`st-ketcher`)
- 🔬 **Descriptor Computation** — Automatically computes molecular descriptors using RDKit
- 🤖 **pIC50 Prediction** — Predicts bioactivity (pIC50) using a pre-trained RASAR (Read-Across Structure–Activity Relationship) model
- 📐 **Applicability Domain (AD) Analysis** — Evaluates whether a query molecule falls within the model's applicability domain using the **leverage approach** (Williams plot)
- 📊 **Result Visualization** — Displays predictions alongside AD warnings, leverage values (*h*), and standardized residuals
- 📥 **Export Results** — Download prediction results as a CSV file

---

## 🖥️ Demo

> Launch the app and draw or paste a SMILES string to get an instant pIC50 prediction with AD analysis.

---

## 🚀 Installation

### Prerequisites

- Python ≥ 3.8
- pip

### Clone the Repository

```bash
git clone https://github.com/yourusername/LipIC-O-Pred.git
cd LipIC-O-Pred
```

### Install Dependencies

```bash
pip install -r requirements.txt
```

### Run the App

```bash
streamlit run app.py
```

The app will open in your browser at `http://localhost:8501`.

---

## 📦 Dependencies

| Package | Role |
|---|---|
| `streamlit` | Web application framework |
| `st-ketcher` | Embedded molecular structure editor (Ketcher) |
| `rdkit` | Molecular descriptor computation & cheminformatics |
| `rasar` | Read-Across Structure–Activity Relationship modelling |
| `pandas` | Data manipulation and result handling |
| `numpy` | Numerical computation (leverage matrix, residuals) |
| `scikit-learn` | Model training/loading utilities |

Install all dependencies via:

```bash
pip install -r requirements.txt
```

A minimal `requirements.txt`:

```
streamlit
streamlit-ketcher
rdkit-pypi
rasar
pandas
numpy
scikit-learn
matplotlib
plotly
```

---

## 🧠 How It Works

### 1. Molecule Input
Users draw a molecule using the Ketcher editor embedded via `st-ketcher`, or paste a SMILES string directly, or give SMILES through uploading the .xlsx file.

### 2. Descriptor Calculation
RDKit computes a set of molecular descriptors that serve as input features for the RASAR descriptor computation.

### 3. RASAR Descriptor Computation
The rasar descriptor (RA function) will be computed based on RDKit descriptors using the Laplacian Kernel similarity function (CTC=4, gamma=0.25). 

### 4. pIC50 Prediction
The pre-trained RASAR model takes the computed descriptors and returns a predicted **pIC50** value.

### 5. Applicability Domain — Leverage Approach
The applicability domain is assessed using the **leverage** (hat value) method:

- The **leverage** *h* of a query compound is computed as:

$$h_i = \mathbf{x}_i^T \left(\mathbf{X}^T \mathbf{X}\right)^{-1} \mathbf{x}_i$$

where **X** is the descriptor matrix of the training set and **xᵢ** is the descriptor vector of the query molecule.

- The **warning leverage** *h** is defined as:

$$h^* = \frac{3(p+1)}{n}$$

where *p* is the number of descriptors and *n* is the number of training compounds.

- A molecule is considered **within the applicability domain** if *h < h**.
- A **Williams plot** (standardized residuals vs. leverage) is generated to visually communicate the AD boundary.

| Condition | Interpretation |
|---|---|
| *h* < *h** | ✅ Within applicability domain |
| *h* ≥ *h** | ⚠️ Outside applicability domain — prediction unreliable |

---

## 📁 Project Structure

```
LipIC-O-Pred/
│
├── app.py                  # Main Streamlit application
├── lib/
│   ├── leverage.py      
│   ├── qsar.xlsx
│   ├── rasar.xlsx       
│   └── logo.png    
├── requirements.txt
├── LICENSE
└── README.md
```

---

## 📊 Output

For each submitted molecule, LipIC-O-Pred returns:

- **Predicted pIC50** value
- **Applicability Domain status** (Inside / Outside)

---

## ⚠️ Limitations

- Predictions are most reliable for molecules structurally similar to the training set.
- Molecules flagged as **outside the applicability domain** (h ≥ h*) should be interpreted with caution.
- The tool is intended for **research and educational purposes** only.

---


## 🤝 Contributing

Contributions, bug reports, and feature requests are welcome! Please open an [issue](https://github.com/yourusername/LipIC-O-Pred/issues) or submit a pull request.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-feature`)
3. Commit your changes (`git commit -m 'Add my feature'`)
4. Push to the branch (`git push origin feature/my-feature`)
5. Open a Pull Request

---

## 📄 License

This project is licensed under the **MIT License** — see the [LICENSE](LICENSE) file for details.

---

## 📬 Contact

For questions or collaborations, please reach out via [GitHub Issues](https://github.com/yourusername/LipIC-O-Pred/issues)

---
