# The following DIR contains:
* 'feature_extraction.py'
  * Process to extract features from DNA and amino acid sequences;
  * It also implements the SurfMap tool - available on (https://github.com/i2bc/SURFMAP?tab=readme-ov-file#Table-of-contents) to
extract surface features of protein structures;


* 'haddock_run.py'
  * Pipeline for subission of aptamer and protein structures in a HADDOCK run;
  * Requires high permission level on the platform - it's free access but requires permission;
 
    
* 'downloader.py'
  * Pipeline that takes a list of links as an input and downloads the best performing structure from
  each HADDOCK run result page;   
