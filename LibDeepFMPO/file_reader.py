"""
This code was proposed by DeepFMPO in:

"Deep reinforcement learning for multiparameter optimization in de novo drug design"
Available at: https://doi.org/10.26434/chemrxiv.7990910.v1

Original repository:
https://github.com/stan-his/DeepFMPO

"""
from LibDeepFMPO.mol_utils import drop_salt
from rdkit import Chem

# Read a file containing SMILES
# The file should be a .smi or a .csv where the first column should contain a SMILES string
def read_file(file_name, drop_first=True):
    
    molObjects = []
    
    with open(file_name) as f:
        for l in f:
            if drop_first:
                drop_first = False
                continue
                
            l = l.strip().split(",")[0]
            smi = drop_salt(l.strip())
            molObjects.append(Chem.MolFromSmiles(smi))
    
    return molObjects

