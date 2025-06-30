# AptaCom Model (V1.5)

This repository contains a Python script for classifying Aptamer-Protein interactions using XGBoost.
This model serves mostly as a preliminary <b>screening tool</b>. The expected usage is to input a CSV file describing 
various aptamer towards one single protein target. Samples are available in the folder 'Sample_runs', feel free to use them as 
reference or even modify them to suit your run. 
Run on python3.8 or higher

## Requirements

To run this script, you'll need the following Python packages:

1.  **XGBoost:** 
    *   Core XGBoost implementation for classification;
    *   Installation: `pip install xgboost`
    *   Installation: `uv add xgboost`

2.  **pandas:**
    *   To load input (CSV shaped file);
    *   Installation: `pip install pandas`
    *   Installation: `uv add pandas`

3.  **pyBioMed:**
    *   Extracts features from sequences used in prediction;
    *   Installation: `pip install pyBioMed`
    *   Installation: `uv add pybiomed`

4.  **rust-sasa-python:**
    *   Extracts features from protein pdb files used in prediction;
    *   Installation: `pip install rust-sasa-python`
    *   Installation: `uv add rust-sasa-python`

5.  **MXfold2:**
    *   Extracts Aptamer secondary structure;
    *   Installation: <href>https://github.com/mxfold/mxfold2</href>
    *   In case of unsucessful installation - alternative methods for secondary structure prediction can be used - it only requires a colum 'SS' on the input CSV file;

## Installation Instructions

1.  Clone this repository: `git clone https://github.com/rpgv/AptaCom.git`
2.  Navigate to the repository directory: `cd AptaCom/Model`
3.  Install the required packages using pip: `pip install -r requirements.txt` or `uv add -r requirements.txt`
    *   Alternatively, install packages individually using `pip install xgboost pandas pyBioMed rust-sasa-python` or `uv add xgboost pandas pyBioMed rust-sasa-python` 
    *   Note: Install rich `pip install rich` or `uv add rich` for colorful outputs;
4.  Create a CSV file following this format:
  <table>
    <thead>
      <tr>
        <th>Aptamer Sequence</th>
        <th>Target Sequence</th>
        <th>PDB_ID</th>
        <th>SS (mxfold2 not installed)</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td>AGUGAUGACUGAC</td>
        <td>MLKAPRKKSMKSKM</td>
        <td>1HAO</td>
        <td>...(.)..</td>
      </tr>
      <tr>
        <td>AGTGATGACTGAC</td>
        <td>MLKAPRKKSMKSKM</td>
        <td>1HAO</td>
        <td>...(())..</td>
      </tr>
    </tbody>
  </table>

</body>
</html>

## Usage

`python AptaCom.py  input_data.csv`

### Input data description

1. If <b>mxfold2</b> is successfully installed:
    *    input_data.csv contains three columns: 'Aptamer Sequence','Target Sequence','PDB_ID'
    *    Note: PDB_ID: describes the name of PDB files with <b>protein atoms only</b>. These files should be in the ~/AptaCom/Model/PDBs directory;
  
2. If <b>mxfold2</b> was not installed successfully: 
    *    input_data.csv contains four columns: 'Aptamer Sequence','Target Sequence','PDB_ID','SS'
    *    Note: PDB_ID: describes the name of PDB files with <b>protein atoms only</b>. These files should be in the ~/AptaCom/Model/PDBs directory;
    *    Note: SS: describes the secondary structure of <b>each aptamer</b> to be classified;

