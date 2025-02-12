import argparse
import csv
import os
from utils import Codon, DFA, NodeType, get_ACGU_char, utr_to_aa, convert_aa_to_triple, prepare_codon_unit_lattice

def get_dfa(aa_graphs, protein, utr_trimmed):
    dfa = DFA()
    newnode = NodeType(3 * len(protein), 0)
    is_utr = -1

    dfa.add_node(newnode)

    for i, aa in enumerate(protein):
        i3 = i * 3
        graph = aa_graphs[aa]

        for pos in range(3):
            for node in graph.nodes[pos]:
                num = node.num
                newnode = NodeType(i3 + pos, num)
                
                dfa.add_node(newnode)
                for edge in graph.right_edges[node]:
                    n2 = edge.node
                    nuc = edge.nuc
                    num = n2.num
                    newn2 = NodeType(i3 + pos + 1, num)
                    if is_utr == -1:
                        dfa.add_edge(newnode, newn2, nuc, round(edge.weight, 3))
                    elif is_utr < len(utr_trimmed):
                        if utr_trimmed[is_utr] == get_ACGU_char(nuc):
                            dfa.add_edge(newnode, newn2, nuc, 0)
                            # dfa.add_edge(newnode, newn2, nuc, round(edge.weight * 100, 3))
                        # else:
                            # dfa.add_edge(newnode, newn2, nuc, 0)
            if is_utr > -1:
                is_utr += 1
        if aa == "STOP" and is_utr == -1:
            is_utr = 0
    return dfa

def dfa_generator(seq, lambda_val=0, utr3=""):
    SEQ = seq
    UTR3 = utr3
    LAMBDA_VAL = lambda_val
    
    current_file_path = os.path.abspath(__file__)
    current_dir = os.path.dirname(current_file_path)
    DFA_FILE = current_dir + "/dfa.txt"
    CODON_TABLE = current_dir[:-8] + "/codon_usage_freq_table_human.csv"
    CODING_WHEEL = current_dir[:-8] + "/coding_wheel.txt"

    codon_table = {}
    with open(CODON_TABLE, mode='r') as codon_file:
        reader = csv.DictReader(codon_file, fieldnames=['codon', 'aa', 'freq'], delimiter=',')
        next(reader)
        for row in reader:
            codon = row['codon']
            aa = row['aa']
            freq = float(row['freq'])
            codon_table[codon] = (aa, freq)

    utr_trimmed = UTR3[:len(UTR3) - (len(UTR3) % 3)]
    utr_aa = utr_to_aa(utr_trimmed, codon_table)

    if UTR3 != "" and SEQ[-1] != '*':
        SEQ = SEQ + '*'
    aa_seq = SEQ + utr_aa

    codon = Codon(CODON_TABLE)
    aa_graphs_with_ln_weights = prepare_codon_unit_lattice(CODING_WHEEL, codon, lambda_=LAMBDA_VAL)
    
    aa_tri_seq = convert_aa_to_triple(aa_seq)
    protein = aa_tri_seq.split()
    dfa = get_dfa(aa_graphs_with_ln_weights, protein, utr_trimmed)

    with open(f'{DFA_FILE}', 'w') as f:
        print(f"{SEQ}", file=f)
        print(f"{utr_trimmed}", file=f)
        print(f"{aa_seq}", file=f)
        dfa.print(file=f)
    
def main():
    parser = argparse.ArgumentParser(description="Generate DFA from sequence")
    parser.add_argument("seq", type=str, help="The amino acid sequence")
    parser.add_argument("-l", "--lambda_val", type=float, default=0, help="The lambda value for calculating edge weight")
    parser.add_argument("-u3", "--utr3", type=str, default="", help="The 3'UTR sequence")

    args = parser.parse_args()
    dfa_generator(args.seq, args.lambda_val, args.utr3)

if __name__ == "__main__":
    main()
