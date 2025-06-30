import os 

path = './PDBs'

files = list(os.listdir(path))

def clean_file(file):
	with open(file, "r") as op:
		data = op.read().split("\n")
	with open(file.replace(".pdb", "_clean.pdb"), "w") as op:
		for i in data:
			if "MASTER" not in i:
				if "ATOM" in i or "TER" in i or "END" in i:
					op.write(f"{i}\n")

if __name__ == "__main__":
	for i in files:
		print(f"Working on {i}")
		clean_file(f"{path}/{i}")

