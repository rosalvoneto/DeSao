import sys
import numpy as np
from LibDeepFMPO.file_reader import read_file
from LibDeepFMPO.mol_utils import get_fragments
from LibDeepFMPO.build_encoding import get_encodings, encode_list, save_decodings
from SimulatedAnnealing import search
from rewards2 import clean_good

import logging
logging.getLogger().setLevel(logging.INFO)


import yaml


def main(fragment_file, lead_file,MaxIterations, temperatura_inicial, alpha, S, s_temp, SampleSize):
    fragment_mols = read_file(fragment_file)
    lead_mols = read_file(lead_file)
    fragment_mols += lead_mols

    logging.info("Read %s molecules for fragmentation library", len(fragment_mols))
    logging.info("Read %s lead moleculs", len(lead_mols))

    fragments, used_mols = get_fragments(fragment_mols)
    logging.info("Num fragments: %s", len(fragments))
    logging.info("Total molecules used: %s", len(used_mols))
    assert len(fragments)
    assert len(used_mols)
    
    print('Get Encodings...')
    encodings, decodings = get_encodings(fragments)
    save_decodings(decodings)
    logging.info("Saved decodings")
    
    lead_mols = np.asarray(fragment_mols[-len(lead_mols):])[used_mols[-len(lead_mols):]]

    X = encode_list(lead_mols, encodings)

    X = clean_good(X, decodings)

    logging.info("Searching...")
    search(X, decodings, MaxIterations, temperatura_inicial, alpha, S, s_temp, SampleSize)

config = yaml.load(open('Config.yaml', 'r'), Loader=yaml.FullLoader)

fragment_file  = config['fragment_file']
print(f'fragment_file: {fragment_file}')

lead_file  = config['lead_file']
print(f'lead_file: {lead_file}')

MaxIterations = config['MaxIterations']
temperatura_inicial = config['temperatura_inicial']
alpha = config['alpha']
S = config['S']
s_temp = config['s_temp']
SampleSize = config['SampleSize']

main(fragment_file, lead_file,MaxIterations, temperatura_inicial, alpha, S, s_temp, SampleSize)