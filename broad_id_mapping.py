#!/usr/bin/env python

from Bio import Entrez
import argparse
import sys

def get_pubchem_cid_dict(mapping_file_path):
    d = {}
    with open(mapping_file_path,mode='rb') as f:
        for line in f.readlines():
            broad_id, cid = line.strip('\n').split('\t')
            d[broad_id] = cid
    return d

def get_smiles(broad_id_path, pubchem_dict, outfile):
    counter = 0
    with open(broad_id_path,mode='rb') as f:
        for line in f.readlines():
            broad_id = line.strip('\n')
            cid = pubchem_dict.get(broad_id,'') # default is empty sring
            if cid == '':
                smiles = ''
            else:
                record = Entrez.read(Entrez.esummary(db="pccompound", id=cid, retmode="xml"))
                if record is not None and len(record) > 0:
                    smiles = record[0].get('CanonicalSmiles','') # default is empty sring
            outfile.write('\t'.join([broad_id, cid, smiles]) + '\n')
            counter += 1
            if counter % 100 == 0:
                sys.stderr.write("Completed %s records.\r"%(counter))
                sys.stdout.flush()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process individuals with variants from a VCF.')
    parser.add_argument('--broad_id_path', '-b', dest='broad_id_path', action='store',
                    help='Path to 1-column text file of Broad IDs')
    parser.add_argument('--mapping_file_path', '-m', dest='mapping_file_path', action='store',
                    help='Path to 2-column text file mapping Broad ID to PubChem ID, from Identifier Exchange Service')
    parser.add_argument('--entrez_email', '-e', dest='entrez_email', action='store', default=None,
                    help='Email address for Entrez API')
    args = parser.parse_args()
    if args.entrez_email is not None:
        Entrez.email = args.entrez_email
    pubchem_dict = get_pubchem_cid_dict(args.mapping_file_path)
    get_smiles(args.broad_id_path, pubchem_dict, sys.stdout)

