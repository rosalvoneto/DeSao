"""
This code was proposed by DeepFMPO in:

"Deep reinforcement learning for multiparameter optimization in de novo drug design"
Available at: https://doi.org/10.26434/chemrxiv.7990910.v1

Original repository:
https://github.com/stan-his/DeepFMPO

We adopted this code without modifications for use in our system.
"""
# GLOBAL PARAMETERS

# Fragmenting and building the encoding
MOL_SPLIT_START = 70
MAX_ATOMS = 12
MAX_FREE = 3
MAX_FRAGMENTS = 12


# Similarity parameters
ETA = 0.1

# Generation parameters
MAX_SWAP = 5
FEATURES = 4


TIMES = 8
