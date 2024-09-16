# In silico aptamer development for potential marine bioremediation

``` Aptamers are single stranded nucleotides, with applications comparable to those of antibodies - targeted drug delivery, bio-sensors etc - while being more chemically malleable.```

In this repo, we compile the tools and data utilized during the development of this thesis. 

* /raw_data/AptaCom -> Consists of the compiled non redundant data from ana array of various aptamer databases:
  * The dataset describes aptamers, their targets and their interaction as follows: 
<table>
  <tr>
    <th>Aptamer (A)</th>
    <th>A.Chemistry</th>
    <th>A.Sequence</th>
    <th>Target (T)</th>
    <th>T.Chemistry</th>
    <th>T.Sequence</th>
    <th>T.External ID</th>
    <th>Affinity</th>
    <th>Reference</th>
    <th>DOI</th>
    <th>Origin</th>
  </tr>
  <tr>
    <td>Name</td>
    <td>DNA/RNA</td>
    <td>'AGCGATGCATG'</td>
    <td>Name</td>
    <td>Protein/Small Molecule</td>
    <td>'MAEVLTLAT'</td>
    <td>PDB/PubChem/ATCC</td>
    <td>0.9pM</td>
    <td>Li.et Al</td>
    <td>----</td>
    <td>UTexas/AptaDB...</td>
  </tr>
</table>

* Model_Predict.json describes a XGBoost model built for predicting aptamer target interactions (Binding vs No Binding)
   * This model is implemented as such:

     ```repo not complete feature model still not up```
     
* Feature_extraction.py describes the pipeline to retrieve aptamer/target (protein) features based on ```full pipeline still not up, functions are working tho```:
      <table>
        <tr>
          <th>Aptamer Sequence</th>
          <th>Target Sequence </th>
          <th>Target PDB File</th>
        </tr>
      </table>

  
