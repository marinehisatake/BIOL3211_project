"""
Author: Hisatake Ishida
Date : 4/5/21
Description: A script that runs multiple sequence alignment on sequence obtained based on imported list of accession
number from excel file.

Incorporates Modules coded by Associate Professor Mikael Boden at the University of Queensland (Email:m.boden@uq.edu.au):
sym: Alphabet setting.
prob: Probabilities setting. (depending on sym)

Incorporates Classes coded by Associate Professor Mikael Boden at the University of Queensland (Email:m.boden@uq.edu.au):
Sequence: Defines a sequence.
Alignment: Defines a multiple sequence alignment.

Incorporates Methods designed by Associate Professor Mikael Boden at the University of Queensland (Email:m.boden@uq.edu.au) to load and run relative files (FASTA),
and to retrieve data from web services (GenBank and ClustalW 2 alignment)
"""

# Importing modules.
import pandas as pd
from sequence import *
import xlrd

# Read Excel file to import list of the accession number of species of interest.
excel_accession_data = pd.read_excel(r'')
philinoidea_id_list = excel_accession_data['Accession'].tolist()

# Creating empty list to extract sequence data from philinoidea_id_list.
seqs_list = []
seq_name_list = []
seq_only_list = []
raw_seq_list = []

# Iterate over lists of accession number.
for i in range(len(philinoidea_id_list)):
    # Storing raw data and creating list of sequence data retrieved based on the accession number from genbank.
    seq_raw = getSequence(philinoidea_id_list[i], 'genbank')
    raw_seq_list.append(seq_raw)
    seq_only_list.append(str(seq_raw.sequence))
    # Extracting the species information (removing unnecessary information).
    old_name = str(seq_raw.name)
    old_info = str(seq_raw.info)
    new_name_with = old_info.replace((old_name+' '), '')
    # Creating a list of species information.
    if ' ' in str(new_name_with):
        words = str(new_name_with).split()
        new_name_seq = words[0]+'_'+words[1]
        loop_count = 0
        while new_name_seq in seq_name_list:
            loop_count = loop_count + 1
            if loop_count == 1:
                new_name_seq = new_name_seq+"_"+str((seq_name_list.count(new_name_seq)+1))
            elif loop_count > 1:
                new_name_seq = new_name_seq[:-1]+str(loop_count)
        if new_name_seq not in seq_name_list:
            seq_name_list.append(new_name_seq)

# Iterating over the species name list to create FASTA file with species information
# and the relevant unaligned sequences.
for index, new_name in enumerate(seq_name_list):
    seq_only_with_name = Sequence(sequence=seq_only_list[index], name=new_name)
    seqs_list.append(seq_only_with_name)
    print(seq_only_with_name)
    print(seq_only_with_name.alphabet)

# Running ClustalW 2 multiple sequence alignment.
multiple_sequence_alignment_clustal = runClustal(seqs_list)

# Writing multiple sequence alignment FASTA file.
writeFastaFile('unaligned_seq.fa', seqs_list)
writeFastaFile('raw_seq.fa', raw_seq_list)
writeFastaFile('aligned_seq.fa', multiple_sequence_alignment_clustal)