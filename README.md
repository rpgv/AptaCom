<h1><img width="50" height="50" alt="Logo Whiteboard" src="https://github.com/user-attachments/assets/d164d0be-eff1-4095-a6f3-57256b6500ba" />AptaCom</h1>


<h3>In silico aptamer development for potential marine bioremediation</h3>
<p>Aptamers are single stranded nucleotides, with applications comparable to those of antibodies - targeted drug delivery, bio-sensors etc - while being more chemically malleable.</p>

<p>In this repo, we compile the tools and data utilized during the development of the work: 'Designing in silico aptamers for potential marine bioremediation'.</p>

* ```/Database/AptaCom``` Contains a non redundant array of data from various aptamer databases:
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

* ```/Model``` describes a XGBoost model built for predicting aptamer target interactions (Binding vs No Binding)
   * For more details: <href>https://github.com/rpgv/AptaCom/tree/main/Model</href>



  
