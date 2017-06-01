import itertools


def reverse_complement(sequence):
    rev_comp = ""
    Watson_Crick = {"A": "T", "C": "G", "T": "A", "G": "C", "a": "t", "t": "a", "c": "g", "g": "c"}
    for base in sequence:
        rev_comp = Watson_Crick[base] + rev_comp
    return rev_comp


# optimizes the overhangs used in Golden Gate assembly
def golden_gate_optimization(sequence_list, backbone_sequence):
    seq_pairs = {0: "ccct", 1: "gctc", 2: "cggt", 3: "gtgc", 4: "agcg", 5: "ctgt", 6: "tgct", 7: "atgg", 8: "gact",
                 9: "ggac", 10: "tccg", 11: "ccag", 12: "cagc", 13: "gttg", 14: "cgaa", 15: "ccat"}
    seq_matches = []
    for x in range(len(sequence_list)+1):
        seq_matches.append([])
        if x == 0:
            for y in range(16):
                if seq_pairs[y] in reverse_complement(backbone_sequence[-35:]):
                    seq_matches[x].append(seq_pairs[y])
                elif seq_pairs[y] in sequence_list[x][:35]:
                    seq_matches[x].append(seq_pairs[y])
        elif x == len(sequence_list):
            for y in range(16):
                if seq_pairs[y] in reverse_complement(sequence_list[x-1][-35:]):
                    seq_matches[x].append(seq_pairs[y])
                elif seq_pairs[y] in backbone_sequence[:35]:
                    seq_matches[x].append(seq_pairs[y])
        else:
            for y in range(16):
                if seq_pairs[y] in reverse_complement(sequence_list[x][-35:]):
                    seq_matches[x].append(seq_pairs[y])
                elif seq_pairs[y] in sequence_list[x][:35]:
                    seq_matches[x].append(seq_pairs[y])
    combs = []
    for x in itertools.product(*seq_matches):
        combs.append(x)
    for comb in combs:
        if len(comb) == len(set(comb)) and len(comb) == len(sequence_list)+1:
            return comb
    # if there are no possible combinations
    return None


#print(golden_gate_optimization(["ccctgtgcccctgtgcccctgtgccagggcctgtgc", "gtgcagcggtgcagcggtgcagcggtgcagcggtgcagggagcg",
#                                "atggatggatggatggatggatggatggatggatggatggatggatggagggatgg",
#                                "gactccatgactgactgactgactgactgactgactgactgactgactaggggact"],
#                               "ggggccatcgaagggggggggggggggggggggggggggggggggggggggggggggggggggggggatggaggggg"))

golden_gate_overhangs = [
                "ccct","gctc","cggt","gtgc","agcg","ctgt","tgct","atgg","gact","ggac","tccg","ccag","cagc","gttg","cgaa","ccat"
            ]

'''
if golden_gate_method == "regular_assembly":
    for unpacked_list in builds_list:
        for part in unpacked_list:
            part.assembly_method = "Type_II_Restriction_Enzyme"
            if "ggtctc" or "gagacc" in part.sequence.lower():
                golden_gate_error = "BsaI_in_seq"
            for i in range(len(unpacked_list)):
                if golden_gate_error != "":
                    break
                if len(unpacked_list) > 16:
                    break
                unpacked_list[i].primer_forward = "aaggtctca" + golden_gate_overhangs[i] + unpacked_list[i].sequence[:20]
                unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-20:] + golden_gate_overhangs[i+1] + "agagaccaa")
                unpacked_list[i].sequence = "aaggtctca" + golden_gate_overhangs[i] + unpacked_list[i].sequence + reverse_complement(golden_gate_overhangs[i+1]) + "agagaccaa"
        new_backbone_sequence = "aaggtctca" + golden_gate_overhangs[len(unpacked_list)+1] + backbone_sequence + reverse_complement(golden_gate_overhangs[0]) + agagaccaa"
if golden_gate_method == "scarless_assembly":
    for unpacked_list in builds_list:
        for part in unpacked_list:
            part.assembly_method = "Type_II_Restriction_Enzyme"
            if "ggtctc" or "gagacc" in part.sequence.lower():
                golden_gate_error = "BsaI_in_seq"
        new_backbone_sequence = backbone_sequence
        for i in range(len(unpacked_list)):
            if golden_gate_error != "":
                break
            if len(unpacked_list) > 16:
                break
            gg_opt = golden_gate_optimization(unpacked_list,backbone_sequence)
            if gg_opt is None:
                golden_gate_error = "no_efficient_overhang_combinations"
                break
            if i == 0:
                if gg_opt[0] in reverse_complement(backbone_sequence[-35:]):
                    new_backbone_sequence = backbone_sequence[:backbone_sequence.find(reverse_complement(gg_opt[0]))] + reverse_complement(gg_opt[0]) + "agagaccaa"
                    unpacked_list[0].primer_forward = "aaggtctca" + backbone_sequence[backbone_sequence.find(reverse_complement(gg_opt[0])):] + unpacked_list[0].sequence[:20]
                    unpacked_list[0].sequence = "aaggtctca" + backbone_sequence[backbone_sequence.find(reverse_complement(gg_opt[0])):] + unpacked_list[0].sequence
                if gg_opt[0] in unpacked_list[0].sequence[:35]:
                    new_backbone_sequence = backbone_sequence + unpacked_list[0].sequence[:unpacked_list[0].sequence.find(gg_opt[0])] + reverse_complement(gg_opt[0]) + "agagaccaa"
                    unpacked_list[0].primer_forward = "aaggtctca" + gg_opt[0] + unpacked_list[0].sequence[unpacked_list[0].sequence.find(gg_opt[0]):unpacked_list[0].sequence.find(gg_opt[0])+20]
                    unpacked_list[0].sequence = "aaggtctca" + unpacked_list[0].sequence[unpacked_list[0].sequence.find(gg_opt[0]):]
'''

import xml.etree.ElementTree as ET
import io
with open()
tree = ET.parse(io.StringIO(xml))
root = tree.getroot()
print(root)