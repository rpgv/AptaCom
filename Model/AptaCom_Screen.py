#
#  CLI Interface for model usage and testing - for now only mac compatible;
# 
from PyBioMed.PyDNA.PyDNAac import *
from PyBioMed.PyDNA.PyDNAac import GetDAC
from PyBioMed.PyDNA import PyDNApsenac
from PyBioMed.PyDNA.PyDNApsenac import GetPseDNC
from PyBioMed.PyDNA.PyDNAnac import *
from PyBioMed.PyDNA.PyDNApsenac import *
from PyBioMed.PyDNA.PyDNAutil import *
from PyBioMed import Pyprotein
import pandas as pd 
from rich import print 
from rust_sasa_python import calculate_sasa_at_residue_level
import subprocess
import xgboost
import numpy as np
import sys


#
# First - Pipeline from ExtractFeature.py - To retrieve input features
#
############# PyBioMed Feature Extraction #######################

def GetAll_Pair(dna,prot):
    resdna = {}
    resdna.update(GetDAC(dna))
    resdna.update(GetDCC(dna))
    resdna.update(GetDACC(dna))
    resdna.update(GetTAC(dna))
    resdna.update(GetTACC(dna))
    resdna.update(GetTCC(dna))
    resdna.update(GetKmer(dna))
    resdna.update(GetPseDNC(dna))
    resdna.update(GetPseKNC(dna))
    resdna.update(GetSCPseDNC(dna))
    resdna.update(GetSCPseTNC(dna))

    resprot = Pyprotein.PyProtein(prot).GetALL()
    return resdna, resprot


def empty_dicts():
    # The resulting dicts from the GetAll_Pair function doesn't match that of a DFrame;
    # To convert it I'll create the following DICT pair - with reference random DNA/Protein set;
    ex_apt = "AGTCGATGGCTGAGGGATCGATG"
    ex_trg = "MWLGRRALCALVLLLACASLGLLYASTRDAPGLRLPLAPWAPPQSPRRVTLTGEGQADLTLVSLDESQMAKHRLLFFKHRLQCMTSQ"
    DNA_Frame, PROT_Frame = GetAll_Pair(ex_apt, ex_trg)
    dna_f = {}
    prot_f = {}
    for p in DNA_Frame.keys():
        # DNA_Frame[f"apt_{p}"] = []
        dna_f[f"apt_{p}"] = []
    for p in PROT_Frame.keys():
        # PROT_Frame[f"trgt_{p}"] = []
        prot_f[f"trgt_{p}"] = []
    return dna_f, prot_f

def get_features(df):
    dna_frame, prot_frame = empty_dicts()
    dna_frame["Aptamer Sequence"] = []
    prot_frame["Aptamer Sequence"] = []
    dna_frame["Target Sequence"] = []
    prot_frame["Target Sequence"] = []
    dna_frame["PDB_ID"] = []
    prot_frame["PDB_ID"] = []
    for n in range(len(df)):
        print(f"Iteration: {n}/{len(df)}")
        print(str(df["Aptamer Sequence"].iloc[n]))
        print(str(df["Target Sequence"].iloc[n]))
        apt = str(df["Aptamer Sequence"].iloc[n]).replace("U", "T").replace("-", "").replace(".", "").replace("3", "")
        og = str(df["Aptamer Sequence"].iloc[n])
        trgt = str(df["Target Sequence"].iloc[n])
        pdb = str(df["PDB_ID"].iloc[n])
        dna_frame["Aptamer Sequence"].append(apt)
        prot_frame["Aptamer Sequence"].append(og)
        dna_frame["Target Sequence"].append(trgt)
        prot_frame["Target Sequence"].append(trgt)
        dna_frame["PDB_ID"].append(pdb)
        prot_frame["PDB_ID"].append(pdb)
        pdna, pprot = GetAll_Pair(apt, trgt)
        for k, v in pdna.items():
            dna_frame[f"apt_{k}"].append(v)
        for k, v in pprot.items():
            prot_frame[f"trgt_{k}"].append(v)
    dna_df = pd.DataFrame.from_dict(dna_frame)
    prot_df = pd.DataFrame.from_dict(prot_frame)
    return dna_df, prot_df


#################################################################


############ Protein SASA Extraction ############################
aa_list = ["GLY","ALA","VAL","LEU","ILE","THR","SER","MET","CYS","PRO","PHE","TYR","TRP","HIS","LYS","ARG","ASP","GLU","ASN","GLN"]
def residue_exposure_map(pdb_file_path): #calculates sasa for each restype in protein; 
    sasa = calculate_sasa_at_residue_level(pdb_file_path) 
    print(sasa)
    sasa_per_restype = {f"sasa_{i}":0 for i in aa_list} # replace the list of nn by list of aa
    for i in sasa:
        restype = str(f"sasa_{i[0].split('_')[1]}")
        sasa_per_restype[restype]+=float(i[1])
    return sasa_per_restype

def build_df(dataframe:pd.DataFrame):#This function receives a dataframe containing a aptamer/target features and identification keys relating to pdb files;
    # Also receives a path where the files are stored
    # tid = [dataframe["PDB_ID"].iloc[i][1] for i in range(len(dataframe))]
    tid = list(dataframe["PDB_ID"])
    print(tid)
    sasa_df = {f"sasa_{i}":[] for i in aa_list} # replace the list of nn by the list of aa
    # aptamer_id = df["Entry_ID"]
    l = len(tid)
    c = 0
    for i in tid[1:]:
        print(f"Calculating SASA for {i}.pdb | Index: {c} of {l}")
        c += 1
        sasa_per_restype =residue_exposure_map(f"./PDBs/{i}_structure_clean.pdb") # Here pdb sub-dir path can be modified
        for j in list(sasa_per_restype.keys()):
            sasa_df[j].append(sasa_per_restype[j])
    return pd.DataFrame().from_dict(sasa_df)


#################################################################


############ Aptamer SS Extraction ##############################

def predict_ss(bs):
    ss_table = {'Aptamer Sequence':[], 'SS':[]}
    for n in range(len(bs)):
        with open('test.fasta', 'w') as op:
            op.write(f">Aptamer\n{bs['Aptamer Sequence'][n]}")    
        out = subprocess.run(["mxfold2", "predict", "test.fasta"], capture_output=True).stdout.decode("utf-8")
        print(n, out.split("\n")[2].split(" ")[0], "- [.]:",str(out.split("\n")[2].split(" ")[0]).count("."), "- [(]:",str(out.split("\n")[2].split(" ")[0]).count("("), "- [)]:",str(out.split("\n")[2].split(" ")[0]).count(")"))
        ss_table["Aptamer Sequence"].append(bs['Aptamer Sequence'][n])
        ss_table["SS"].append(out.split("\n")[2].split(" ")[0])
    # if bs["Aptamer Sequence"][n] == stop_seq:
    #     break
    return ss_table


def analyse_ss(df):
    ss_table = {"Aptamer Sequence":df["Aptamer Sequence"], "SS":df["SS"]}
    ss_table["SG_."] = []
    ss_table["SG_("] = []
    ss_table["SG_)"] = []
    ss_table["S1_."] = []
    ss_table["S1_("] = []
    ss_table["S1_)"] = []
    ss_table["S2_."] = []
    ss_table["S2_("] = []
    ss_table["S2_)"] = []
    ss_table["S3_."] = []
    ss_table["S3_("] = []
    ss_table["S3_)"] = []
    ss_table["S4_."] = []
    ss_table["S4_("] = []
    ss_table["S4_)"] = []
    sc = 0
    for ss in ss_table["SS"]:
        l = len(ss)
        d = len(ss)//4
        r = len(ss)%4
        c = 0
        
        ss_table["SG_."].append(ss.count("."))
        ss_table["SG_("].append(ss.count("("))
        ss_table["SG_)"].append(ss.count(")"))   
    
        for s, i in zip(["S1","S2","S3", "S4"],range(0,l-d,d-1)):
            c+=1
            print("----------")
            ss_table[f"{s}_."].append(ss[i:i+d].count("."))
            print(i,i+d,d)
            ss_table[f"{s}_("].append(ss[i:i+d].count("("))
            ss_table[f"{s}_)"].append(ss[i:i+d].count(")"))
            print("----------")
        if c != 4:
            sc += 1
            print(f"SS: {ss} | Len: {len(ss)} | ")
    
        c = 0
    print(pd.DataFrame().from_dict(ss_table).shape, pd.DataFrame().from_dict(ss_table).columns)
    ndf = pd.DataFrame().from_dict(ss_table)
    return ndf

#################################################################

if __name__ == "__main__":
    print("[blue]Starting Feature Extraction Pipeline[/blue]")
    try:
        df_file = sys.argv[1] # CSV file with Aptamer and Protein Sequence and Protein PDB ID
    except: 
        print("[bold red]Missing arguments[/bold red]")
        print("[green]Correct usage: ExtractFeatures.py <CSV.file>[/green]\nCSV.file describes:\n-Aptamer sequence\n-Protein sequence\n-Protein: PDB_ID")
        print("Script assumes protein structures are in sub-directory ./PDBs\n")

    ############## Loading CSV #################################
    print("Loading CSV file...")
    data = pd.read_csv(df_file)
    print("File loaded successfully!")
    print(f"File - {df_file} - of shape: {data.shape}")
    print("\n")
    ############### Starting PyBioMed pipeline ##################
    print("Starting PyBioMed Token and Physiochemical features")
    data["Target Sequence"] = data["Target Sequence"].apply(lambda x: x[2:-2].split(',')) # Some sequences come with formatting issues i.e: ['MWRTLP']
    #Due to some formatting issues, there are some corrections to be made to sequences;
    target_sequence = []
    for i in range(len(data["Target Sequence"])):
        target_sequence.append(str(list(data["Target Sequence"].iloc[[i]])[0][0]))
    data["Target Sequence"] = target_sequence
    data["Aptamer Sequence"] = data["Aptamer Sequence"]
    data.dropna(inplace=True)
    print(f"Data Shape after dropna and other modifications: {data.shape}")
    print("Begining calculations...")
    dna_df, prot_df = get_features(data) #Here we extract the features of both protein and aptamer; 
    print(f"Checking output shapes - Aptamer: {dna_df.shape} and Target: {prot_df.shape}")
    merged_token = pd.concat([dna_df,prot_df], axis=1)
    print(f"Checking Shapes after merging: Feature: {merged_token.shape}")
    merged_token.to_csv("TokenFeatures_screening.csv", index= False) # Intermediate file
    print("[green]FINISHED PYBIOMED[/green]")
    ############# Finished PyBioMed pipeline ####################


    ############ Starting Protein SASA Extraction ###############
    print("Starting with proein SASA calculation (per residue)")
    merged_token = pd.read_csv("TokenFeatures_screening.csv")
    output = build_df(merged_token)
    print(f"Protein SASA Shape: {output.shape}")
    merged_sasa = pd.concat([merged_token, output], axis=1)
    merged_sasa.to_csv("ProteinSasa_screening.csv",index=False) # Intermediate file
    print("[green]FINISHED SASA calculation[/green]")
    ########### Finished Protein SASA Extraction  ###############


    ########## Starting Aptamer SS prediction ###################
    try:
        print("Starting aptamer Secondary Structure prediction")
        merged_sasa = pd.read_csv("ProteinSasa_screening.csv")
        ss_table = predict_ss(merged_sasa)
        ss_df = pd.DataFrame().from_dict(ss_table)
        print(f"Aptamer SS_df table shape: {ss_df.shape}")
        print("Segmenting SS counts per aptamer section - divides apt sequence into 4 segments and accounts for number of ss units per section")
        #ss_df.to_csv("Aptamer_Secondary_Structure_Table.csv", index=False) # Intermediate file
        ss_segments = analyse_ss(ss_df)
        print("Segmented SS counts shape: {ss_segments.shape}")
        merged_ss = pd.concat([ss_segments, merged_sasa], axis=1)
        ########## Finished Aptamer SS prediction ###################
        print("[green]Process finished saving file to ExtractedFeatures.csv ...[/green]")
        merged_ss.to_csv("ExtractedFeatures_screening.csv", index=False)
    except Exception as e:
        print("[red]\nError in Aptamer Secondary Structure Prediction\nMight be due to problems in mxfold2 installation[/red]")
        print(e)
        if 'SS' in list(data.columns):
            print("Falling back on SS structures present in CSV input file")
            ss_segments = analyse_ss(ss_df)
            print("Segmented SS counts shape: {ss_segments.shape}")
            merged_ss = pd.concat([ss_segments, merged_sasa], axis=1)
            print("[green]Process finished saving file to ExtractedFeatures.csv ...[/green]")
            merged_ss.to_csv("ExtractedFeatures_screening.csv", index=False)
        else:
            print("[red]Could not complete feature extraction due to missing Secondary Strucutre column and missing mxfold2 installation[/red]")
            quit()

    ############## Concatenating Datasets #######################

    

    #############################################################


    ############# Starting Model Prediction #####################
    # (Load Model)
    m3 = xgboost.XGBClassifier()
    m3.load_model("Model_Class_m16.json")
    # Load Screening dataset
    df = pd.read_csv("./ExtractedFeatures_screening.csv")
    y_true = df["Binding"]
    to_drop = [i for i in df.columns if "Aptamer" in i or "Target" in i or "Entry" in i or "Binding" in i or i == "SS" or "PDB_ID" in i]
    to_drop = to_drop + ["Binding"]
    X = df.drop(to_drop, axis=1)
    y_pred = m3.predict(X)
    print(f"Aptamers expected to interact with protein target:")
    binding_aptamers = {"Aptamer Sequence":[]}
    for i, j in zip(df["Aptamer Sequence"], y_pred):
        if j == 1:
            print(f"Aptamer: {i} predicted to [green]bind[/green]")
            binding_aptamers["Aptamer Sequence"].append(i)
    out = pd.DataFrame().from_dict(binding_aptamers)
    out.to_csv("PredictedAptamers.csv", index=False)
    print("Successfully saved predicted aptamers to PredictedAptamers.csv")
        

