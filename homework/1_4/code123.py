from collections import defaultdict
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction

import numpy as np
import random
# from Bio.Seq import Seq
islands_p = "data/islands.fasta"
genome_p = "data/nonIslands.fasta"

def decode_pair_percentage(dictionary='keep'):
    if dictionary == 'keep':
        dictionary = {'AA': 187137, 'AC': 167700, 'AG': 262824, 'AT': 122114, 'CA': 229823, 'CC': 375624, 'CG': 215537, 'CT': 262183, 'GA': 218792, 'GC': 322299, 'GG': 378421, 'GT': 170499, 'TA': 103975, 'TC': 217531, 'TG': 233253, 'TT': 186460}
    if dictionary == 'genome':
        dictionary = {'AA': 19089, 'AC': 10010, 'AG': 14908, 'AT': 17437, 'CA': 14482, 'CC': 10623, 'CG': 1121, 'CT': 15432, 'GA': 12793, 'GC': 7270, 'GG': 10644, 'GT': 11442, 'TA': 15081, 'TC': 13754, 'TG': 15476, 'TT': 22830}

    # print(dictionary)
    
    total = 0
    for key in dictionary:
        total += dictionary[key]
    for key in dictionary:
        dictionary[key] = round(dictionary[key] / total, 4)
    
    let_to_n = {'A': 0, 'C':1, 'G':2, 'T':3}

    dicty = defaultdict(lambda: [0, 0, 0, 0])
    for pair, cnt in dictionary.items():
        fst = pair[0]
        snd = pair[1]
        dicty[fst][let_to_n[snd]] = cnt
    
    # ValueError: probabilities do not sum to 1     ok...
    for letter, vals in dicty.items():
        sum = np.sum(vals)
        for i in range(4):
            dicty[letter][i] = vals[i] / sum

    # print(dicty['A'][1], dictionary['AC'])


    return dicty


def percentage(path):

    nucleotide_cnt={'A': 0, 'C': 0, 'G': 0, 'T': 0}
    total = 0

    fasta_sequences = SeqIO.parse(path, "fasta")
    for record in fasta_sequences: #record's a string
        total+=len(record)
        for nuc in record:
            nucleotide_cnt[nuc] +=1

    print("\nscanning", path)
    for nuc, cnt in nucleotide_cnt.items():
        print(f"{nuc} was found in {round(cnt/total*100)}% nucleotides")
    # print("--finished calculations--\n")

def pair_percentage(path, verbose = 0):

    pair_cnt=dict()
    for n1 in "ACGT":
        for n2 in "ACGT":
            pair_cnt[n1+n2]=0
    total = 0

    fasta_sequences = SeqIO.parse(path, "fasta")
    for record in fasta_sequences: #record's a string
        total+=len(record)-1
        for i in range(1,len(record)):
            pair = record[i-1] + record[i]
            pair_cnt[pair] +=1
    
    if verbose:
        print("\nscanning", path)
        for pair, cnt in pair_cnt.items():
            print(f"pair {pair} was found in rounded {round(cnt/total*100)}% pairs")
        
        print(f"precisely, percent of GC pairs was {round(pair_cnt["GC"]/total*100, 1)}%, \ncompared to expected {round(100/16,1)}%.")
        print("--finished scanning file--\n")
    
    return decode_pair_percentage(pair_cnt)







# islands_pairs = pair_percentage(islands_p)


def hidden_markov(island_probas, genome_probas, n=10, out='genseqs.txt'): ##generates a sequence of len n.

    #chances of States
    isle2isle = 0.95
    isle2gen = 1 - isle2isle
    gen2gen = 0.995
    gen2island = 1-gen2gen


    # island_probas = decode_pair_percentage('keep')
    # genome_probas = decode_pair_percentage('genome')
    # print(island_probas['G']['C'])
    # print(genome_probas['G']['C'])

    
    #random choice of init state

    state = np.random.choice([1, 0])
    choices = ['A', 'C', 'G', 'T']
    proba_first = [0.2,0.3,0.3,0.2] if state == 1 else [0.3,0.2,0.2,0.3]
    
    letter = np.random.choice(choices, p=proba_first)
    Sequence = [letter] #first letter
    # print(Sequence)
    States = [state]

    for i in range(1,n): #first letter
        if state == 1:
            probas = island_probas
            state = np.random.choice([1, 0], p=[isle2isle, isle2gen]) #stay or switch
        else:#genome zero 0
            probas = genome_probas
            state = np.random.choice([0, 1], p=[gen2gen, gen2island]) #stay or switch

        letter = np.random.choice(choices, p=probas[letter]) # dict with lists corresponding to ACGT
        Sequence.append(letter)
        States.append(state)
    
    normal_states = [int(x) for x in States]
    sequence_str = ''.join(Sequence)  


    # return sequence_str, normal_states

    seq_record = SeqRecord(
        Seq(''.join(Sequence)),
        id=f"Generated_sequence_of_len_{n}", description=f'gc_fraction = {gc_fraction(Sequence)}')


    with open(out, 'a') as f:
        SeqIO.write(seq_record, f, "fasta")
        f.write(f"{''.join(map(str, States))}\n\n")  # Add states and blank line
    print(f"Generated sequence of len {n}.")
    

    #     #select random letter. decide if to switch probabilities
    #     #if decided to switch, update probabilities













# from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord

# # Create a SeqRecord with a specific identifier and description
# record = SeqRecord(Seq("ATCGTAGCTAG"), id="sequence_1", description="Example DNA sequence")
# print(record.description)
# # Write the SeqRecord to a FASTA file
# SeqIO.write(record, "example.fasta", "fasta")
# print("Sequence written to example.fasta")

































def read_and_print_fasta_GENERAL(path, cnt, length):
    fasta_sequences = SeqIO.parse(path, "fasta")
    t = 0
    for record in fasta_sequences:
        if t<cnt:
            print("Read sequence from file:")
            print("ID:", record.id)
            print("length: ", len(record)) 
            print("Sequence clip:", record.seq[0:length])
            t+=1
        else:
            print("----examples finished---")
            break

