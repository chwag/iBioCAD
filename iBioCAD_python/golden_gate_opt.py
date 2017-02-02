import itertools

def reverse_complement(sequence):
    rev_comp = ""
    Watson_Crick = {"A":"T","C":"G","T":"A","G":"C","a":"t","t":"a","c":"g","g":"c"}
    for base in sequence:
        rev_comp = Watson_Crick[base] + rev_comp
    return rev_comp

def golden_gate_optimization(sequence_list):
    seq_pairs = {0:"ccct",1:"gctc",2:"cggt",3:"gtgc",4:"agcg",5:"ctgt",6:"tgct",7:"atgg",8:"gact",9:"ggac",10:"tccg",11:"ccag",12:"cagc",13:"gttg",14:"cgaa",15:"ccat"}
    seq_matches = []
    for x in range(len(sequence_list)-1):
        seq_matches.append([])
        for y in range(16):
            if seq_pairs[y] in sequence_list[x][-35:]:
                    seq_matches[x].append(seq_pairs[y])
            if seq_pairs[y] in reverse_complement(sequence_list[x+1][:35]):
                    seq_matches[x].append("-"+seq_pairs[y])
    combs = []
    for x in itertools.product(*seq_matches):
        combs.append(x)
    for comb in combs:
        if len(comb) == len(set(comb)):
            return comb
    #if there are no possible combinations
    return None


print(golden_gate_optimization(["ccctgtgcccctgtgcccctgtgcccctgtgc","gtgcagcggtgcagcggtgcagcggtgcagcggtgcagcg","atggatggatggatggatggatggatggatggatggatggatggatggatgg","gactgactgactgactgactgactgactgactgactgactgactgact"]))