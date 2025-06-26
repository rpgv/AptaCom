# The following DIR contains:
* 'FeatureExtract.py'
  * Process to extract features from DNA and amino acid sequences;
  * It implements the Rust-SASA tool - available on ([https://github.com/i2bc/SURFMAP?tab=readme-ov-file#Table-of-contents](https://github.com/maxall41/rust-sasa-python)) to
extract SASA information from proteins;
  * It implements mxfold2 for Aptamer Secondary Structure prediction - available for installation on (https://github.com/mxfold/mxfold2)
<br>
[Usage: FeatureExtract.py <Apt_Trgt_pair.csv>] Apt_Trgt_pair.csv describes a dataframe with aptamer sequence, protein sequence, and PDB_ID, which describes the path of the respective protein PDB file
<br>
(for this pipeline the pdb file should be clean of non protein data, including connections, remarks, master, and atom charge information)

  
* 'haddock_run.py'
  * Pipeline for subission of aptamer and protein structures in a HADDOCK run;
  * Requires high permission level on the platform - it's free access but requires permission;
 
    
* 'downloader.py'
  * Pipeline that takes a list of links as an input and downloads the best performing structure from
  each HADDOCK run result page;   
