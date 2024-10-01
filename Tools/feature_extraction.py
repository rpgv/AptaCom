#Complete data extraction pipeline
import pandas as pd
import os
import itertools as it
from PyBioMed.PyDNA.PyDNAac import *
from PyBioMed.PyDNA.PyDNAac import GetDAC
from PyBioMed.PyDNA import PyDNApsenac
from PyBioMed.PyDNA.PyDNApsenac import GetPseDNC
from PyBioMed.PyDNA.PyDNAnac import *
from PyBioMed.PyDNA.PyDNApsenac import *
from PyBioMed.PyDNA.PyDNAutil import *
from PyBioMed import Pyprotein
APT = "AGTCGATGGCTGAGGGATCGATG" #These are sample sequences tocreate a template
TRGT = "MWLGRRALCALVLLLACASLGLLYASTRDAPGLRLPLAPWAPPQSPRRVTLTGEGQADLTLLQCMTSQ"
#Path where SurfMap will work from. 
# ON unix systems -> main creates the dirs following the following dir structures:
# Surfmap will retrieve proteins from ./'surfmap_path'/ProteinStructures - structures should be named: 'structure id'_structure.pdb
# Surfmap stores it's results inside ./'surfmap_path'/SurfmapOut/smoothed_matrices
surfmap_path = f"./Run_FeatureExtract"

#This funciton runs trhough pybiomeds features and characterizes the given protein and nucelotide sequence
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

#This funciton utilizes the sample aptamer and target sequences to create a template dicitonary to add features to
def template_dict():
    dnf, tgf = GetAll_Pair(APT, TRGT)
    template_dict = {}
    prot_df = {i:[] for i in tgf.keys()}
    dna_df = {i:[] for i in dnf.keys()}
    template_dict.update(prot_df)
    template_dict.update(dna_df)
    return template_dict

#This function receives a dictionary containing nucleotide sequences and their respective protein targets:
#  sample = {"Aptamer Sequence":["AAAGTGCTGACTGAC"], "Target Sequence": ["MAEVKLPRTAKBG"]}
#  NOTE: PyBioMed only takes DNA sequences. Both sequences MUST be larger than 10 units, without missing or X characters;
def retrieve_features(d:dict):
    apt_sequences = d["Aptamer Sequence"]
    trg_sequences = d["Target Sequence"]
    if len(apt_sequences) == len(trg_sequences):
        for a, t in zip(apt_sequences,trg_sequences):
            print(list(apt_sequences).index(a),"/",len(apt_sequences))
            dna, trg = GetAll_Pair(a, t)
            [d[i].append(j) for i, j in dna.items()]
            [d[i].append(j) for i, j in trg.items()]
    return d #Returns a dict wth the features

#The folowing funtions take as input a datframe. The dataframe 'df'was obtained from the dictionary generated in the previous functions

#This function takes in a dataframe, with a comlum PDB_ID that refers to a protein ids associated with files with the name 'PDB_ID'_structure.pdb
#This implementation of surfmap requires docker to run;
#This function only calculates the circular variance feature of the protein surface
def surfmap_creator(df:pd.DataFrame):
    ids = list(df["PDB_ID"])
    files = os.listdir(f"{surfmap_path}/SurfmapOut/smoothed_matrices/")#Correct paths
    for i in ids:
        if f'{i}_structure_circular_variance_smoothed_matrix.txt' not in files:#Correct paths
            print("\nCurrently creating map for protein:",i)
            os.system(f"surfmap -pdb {surfmap_path}/ProteinStructures/{i}_structure.pdb -tomap circular_variance -d {surfmap_path}/SurfmapOut/ --docker") #Correct paths
        else: pass
#This function doesnt return anything;

#This function receives the same dataframe fed into 'surfmap_creator' function, and uses it to retrieve the output of the surfmap run
#from this txt output it generates a dataframe of surface features organized per PDB_ID - must match filename
def surfmap_retriever(df):
    # output_SURFMAP_{i}_structure_all_properties
    ids = list(df["PDB_ID"])
    ish = df.shape
    merge = pd.DataFrame()
    for i in ids:
        matrix1 = pd.read_csv(
            f'{surfmap_path}/SurfmapOut/smoothed_matrices/{i}_structure_circular_variance_smoothed_matrix.txt',#Correct paths
            delimiter="\t")
        matrix1 = matrix1.rename(columns={"svalue": "svalue1"})
        matrix1 = matrix1["svalue1"]
        merge = pd.concat([merge, matrix1], axis=1)
    merge = merge.transpose()
    merge = merge.reset_index(drop=True)
    extra_columns = [i for i in list(merge.columns) if int(i) >=111 and int(i)<=777]
    new_columns = {i:f"{str(i).replace('.1', '')}_SMap" for i in extra_columns}
    merge.rename(inplace=True, columns=new_columns)

    df = df.reset_index()
    df = pd.concat([df, merge], axis=1)

    print(f'Resulting DF of shape: {df.shape} from original shape (DF): {ish} and original shape (MERGED): {merge.shape}')
    return df #returns dataframe containing surfmap features


#Simple function to retrieve fasta from rcsb
def get_fasta(idx):
    import requests
    url = f'https://www.rcsb.org/fasta/entry/{idx}/display'
    raw = requests.get(url).content
    sequence = 'Default - No sequence found'
    try:
        sequence = str(raw).split("\\n")[1]
    except Exception as e:
        print(f'Exception {e} occurred on id: {id}')
    print(f'{sequence}')
    return sequence


# Converts three letter amino acids into single letter amino acids
def single_letter(a):
    map = {"GLY":	"G",
    "ALA":	"A",
    "VAL":	"V",
    "LEU":	"L",
    "ILE":	"I",
    "THR":	"T",
    "SER":	"S",
    "MET":	"M",
    "CYS":	"C",
    "PRO":	"P",
    "PHE":	"F",
    "TYR":	"Y",
    "TRP":	"W",
    "HIS":	"H",
    "LYS":	"K",
    "ARG":	"R",
    "ASP":	"D",
    "GLU":	"E",
    "ASN":	"N",
    "GLN":	"Q"}
    return map[a]
# Converts single amino acids into six groups of aminoacids:
# Aliphatic (A) / Uncharged (U) / Positiviely charged (P) / Negatively charged (N) / Aromatic (S) / Other (O)
aa_groups = {
    "G" : "A",
    "A" : "A",
    "V" : "A",
    "L" : "A",
    "I" : "A",
    "P" : "O",
    "F" : "S",
    "W" : "S",
    "M" : "A",
    "S" : "U",
    "T" : "U",
    "C" : "U",
    "Y" : "S", 
    "N" : "U",
    "Q" : "U",
    "D" : "N",
    "E" : "N",
    "K" : "P",
    "R" : "P",
    "H" : "U",
}


def convert_to_group(AA:str):
    return aa_groups[AA]  #Employs converstion of single amino acids to group

#This function generates all the possible triad (three) tokens of aminoacids when organized by group 
def get_group_tokens(): 
    return [''.join(i) for i in it.product("PNUASO", repeat=3)] 


#This function calculates three tokens of a given sequence
def getKmers(sequence, size, step):
    sequence = str(sequence)
    return [sequence[x:x+size] for x in range(len(sequence) - size + step)]

# This function generates kmers regarding the group aminocids. It takes a dictionary as input, which should already
# contain a key "target grouped sequence" - in which the amino acid sequence is converted into a grouped anlagous;
def protein_group_kmers(ds:dict):
    tokens = get_group_tokens()
    for t in tokens:
        ds[f'G_{t}']=[]
    for i in ["P","N","U","A","S","O"]:
        ds[f'{i}_GroupContent'] = []
    if "Target Grouped Sequence" not in list(ds.keys()):
        print("Error doesnt contain target grouped sequence column")
        return "Error - Dictionary doesn't contain 'Target Grouped Sequence'"
    else:
        for s in ds["Target Grouped Sequence"]:
            #print("Sequence:", s)
            if s != '':
                s = str(s)
                seq_kmers3 = getKmers(s, 3, 1)
                # print(seq_kmers3)
                for i in ["P","N","U","A","S","O"]:
                    #print("Groups")
                    #print(round(str(s).count(i) / len(s), 3))
                    ds[f"{i}_GroupContent"].append(round(str(s).count(i) / len(s), 3))
                for p in tokens: 
                    #print("Tokens")
                    ds[f'G_{p}'].append(round((seq_kmers3.count(p)/(len(s))), 3)) 
    return ds

#This funtion takes an amino acid sequence and returns a grouped amino acid sequence
def convert_sequence(seq):
    #print("Old:", seq)
    new_seq = ''
    for i in seq:
        new_seq += convert_to_group(i)
    #print("New:",new_seq)
    return new_seq

if __name__ == '__main__':
    os.system("mkdir ./Run_FeatureExtract/SurfmapOut")
    os.system("mkdir ./Run_FeatureExtract/SurfmapOut/smoothed_matrices")
    df = pd.read_csv("./Run_FeatureExtract/screening.csv") #Insert a CSV containing 3 columns Aptamer (DNA) sequence
    # Protein target (amino acid) sequence and protein target ID -> ["Aptamer Sequence", "Target Sequence", "PDB_ID"]
    print(f"Dataframe Shape:{df.shape}, with columns {df.columns}")
    # Converting pandas df to simple dict; 
    base = {i:[] for i in df.columns}
    dfd = df.to_dict()
    for i in dfd:
            base[i] = list(dfd[i].values())
    # Creates a template dict to populate with the converted dict (base)
    template = template_dict()
    template.update(base)
    # Retrieve the sequence related features;
    features = retrieve_features(template)
    # Grouping amino acid sequences 
    features["Target Grouped Sequence"] = []
    for i in template["Target Sequence"]:
        new_sequence = convert_sequence(i)
        features["Target Grouped Sequence"].append(new_sequence)
    new_features = protein_group_kmers(features)
    # Convert to Dataframe
    feat = pd.DataFrame.from_dict(new_features)
    # SurfMap creation - IT CAN'T RUN ON JUPYTER, only IDE
    # surfmap_creator(feat)
    # After running, retrieve structure information
    df = surfmap_retriever(feat)
    df.to_csv("/Users/rpgv2000/Desktop/Masters/Dissertation/Data_analysis/Model/Model_Implementation/sample_aptamer.csv", index=False)
