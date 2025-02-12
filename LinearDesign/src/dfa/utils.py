from collections import defaultdict, namedtuple
import math
import csv

NodeType = namedtuple("NodeType", ["index", "num"])
NodeNucWType = namedtuple("NodeNucWType", ["node", "nuc", "weight"])
NodeNucNodeType = namedtuple("NodeNucNodeType", ["node1", "nuc", "node2"])

### ----------------------------------------------------------------------------
k_map_3_1 = {
    "Phe": 'F',
    "Leu": 'L',
    "Ser": 'S',
    "Tyr": 'Y',
    "STOP": '*',
    "Cys": 'C',
    "Trp": 'W',
    "Pro": 'P',
    "His": 'H',
    "Gln": 'Q',
    "Arg": 'R',
    "Ile": 'I',
    "Met": 'M',
    "Thr": 'T',
    "Asn": 'N',
    "Lys": 'K',
    "Val": 'V',
    "Asp": 'D',
    "Glu": 'E',
    "Gly": 'G',
    "Ala": 'A'
}

### (class) Codon --------------------------------------------------------------
class Codon:
    def __init__(self, path):
        self.codon_table = {}
        self.aa_table = defaultdict(list)
        self.max_aa_table = {}

        with open(path) as codon_file:
            reader = csv.DictReader(codon_file, fieldnames=['codon', 'aa', 'freq'], delimiter=',')
            next(reader)
            for row in reader:
                codon = row['codon']
                aa = row['aa']
                freq = float(row['freq'])
                self.codon_table[codon] = (aa, freq)
                self.aa_table[aa].append((codon, freq))
                self.max_aa_table[aa] = max(self.max_aa_table.get(aa, 0), freq)
            if len(self.codon_table) != 64:
                raise ValueError("Codon frequency file needs to contain 64 codons!")

    def get_weight(self, aa_tri, codon):
        if aa_tri in k_map_3_1:
            aa = k_map_3_1[aa_tri]
            if aa in self.aa_table:
                codons = self.aa_table[aa]
                for c, weight in codons:
                    if c == codon:
                        return weight
        return 0.0

### (class) Lattice ------------------------------------------------------------
class Lattice:
    def __init__(self):
        self.nodes = defaultdict(list)
        self.left_edges = defaultdict(list)
        self.right_edges = defaultdict(list)

    def add_edge(self, n1, n2, nuc, weight=0.0):
        self.right_edges[n1].append(NodeNucWType(n2, nuc, weight))
        self.left_edges[n2].append(NodeNucWType(n1, nuc, weight))

    def add_node(self, n1):
        pos = n1.index
        self.nodes[pos].append(n1)

def prepare_codon_unit_lattice(wheel_path, codon, lambda_=1.0):
    '''
    Step 1: Read initial weights
    Step 2: Update best weights
    Step 3: Read log weights
    '''
    nodes_with_best_weight = defaultdict(lambda: defaultdict(float))
    edges_with_best_weight = defaultdict(lambda: defaultdict(float))
    
    aa_graphs_with_weights = read_wheel_with_weights(wheel_path, nodes_with_best_weight, edges_with_best_weight, codon)
    update_best_weights(aa_graphs_with_weights, nodes_with_best_weight, edges_with_best_weight)
    aa_graphs_with_ln_weights = read_wheel_with_weights_log(wheel_path, nodes_with_best_weight, edges_with_best_weight, codon, lambda_)

    return aa_graphs_with_ln_weights

def read_wheel_with_weights(filename, nodes_with_best_weight, edges_with_best_weight, codon):
    aa_graphs = {}

    with open(filename) as inFile:
        for line in inFile:
            stuff = line.strip().split('\t')
            aa = stuff[0]
            graph = Lattice()
            graph.add_node(NodeType(0, 0))

            last_first = ''
            i = 0
            for option in stuff[1:]:
                option_splited = option.split(' ')
                first = option_splited[0]
                second = option_splited[1]
                thirds = option_splited[2]
                n2 = NodeType(2, i)
                graph.add_node(n2)
                if first != last_first:
                    n1 = NodeType(1, i)
                    graph.add_node(n1)
                    first_num = get_ACGU_num(first)
                    weight = 0.0
                    if NodeType(0, 0) in nodes_with_best_weight.get(aa, {}):
                        weight = edges_with_best_weight[aa][(NodeType(0, 0), first_num, n1)] / nodes_with_best_weight[aa][NodeType(0, 0)]
                    graph.add_edge(NodeType(0, 0), n1, first_num, weight)
                else:
                    n1 = NodeType(1, i-1)

                last_first = first
                second_num = get_ACGU_num(second)
                weight = 0.0

                if n1 in nodes_with_best_weight.get(aa, {}):
                    weight = edges_with_best_weight[aa][(n1, second_num, n2)] / nodes_with_best_weight[aa][n1]
                graph.add_edge(n1, n2, second_num, weight)

                for third in thirds:
                    three_nums = first + second + third
                    weight = 0.0
                    if n2 in nodes_with_best_weight.get(aa, {}):
                        weight = codon.get_weight(aa, three_nums) / nodes_with_best_weight[aa][n2]
                    else:
                        weight = codon.get_weight(aa, three_nums)
                    graph.add_edge(n2, NodeType(0, 0), get_ACGU_num(third), weight)

                i += 1
            aa_graphs[aa] = graph

    return aa_graphs

def update_best_weights(aa_graphs_with_weights, nodes_with_best_weight, edges_with_best_weight):
    for aa, graph in aa_graphs_with_weights.items():
        for node_at_2 in graph.nodes[2]:
            for node_at_3_nuc_weight in graph.right_edges[node_at_2]:
                node_at_3 = node_at_3_nuc_weight.node
                nuc = node_at_3_nuc_weight.nuc
                weight = node_at_3_nuc_weight.weight
                nodes_with_best_weight[aa][node_at_2] = max(nodes_with_best_weight[aa][node_at_2], weight)
                edges_with_best_weight[aa][(node_at_2, nuc, node_at_3)] = weight

        for node_at_1 in graph.nodes[1]:
            for node_at_2_nuc_weight in graph.right_edges[node_at_1]:
                node_at_2 = node_at_2_nuc_weight.node
                nuc = node_at_2_nuc_weight.nuc
                nodes_with_best_weight[aa][node_at_1] = max(nodes_with_best_weight[aa][node_at_1], nodes_with_best_weight[aa][node_at_2])
                edges_with_best_weight[aa][(node_at_1, nuc, node_at_2)] = nodes_with_best_weight[aa][node_at_2]

        for node_at_0 in graph.nodes[0]:
            for node_at_1_nuc_weight in graph.right_edges[node_at_0]:
                node_at_1 = node_at_1_nuc_weight.node
                nuc = node_at_1_nuc_weight.nuc
                nodes_with_best_weight[aa][node_at_0] = max(nodes_with_best_weight[aa][node_at_0], nodes_with_best_weight[aa][node_at_1])
                edges_with_best_weight[aa][(node_at_0, nuc, node_at_1)] = nodes_with_best_weight[aa][node_at_1]


def read_wheel_with_weights_log(filename, nodes_with_best_weight, edges_with_best_weight, codon, lambda_):
    aa_graphs = {}

    with open(filename) as inFile:
        for line in inFile:
            stuff = line.strip().split('\t')
            aa = stuff[0]
            graph = Lattice()
            graph.add_node(NodeType(0, 0))

            last_first = ''
            i = 0
            for option in stuff[1:]:
                option_splited = option.split(' ')
                first = option_splited[0]
                second = option_splited[1]
                thirds = option_splited[2]
                n2 = NodeType(2, i)
                graph.add_node(n2)
                if first != last_first:
                    n1 = NodeType(1, i)
                    graph.add_node(n1)
                    first_num = get_ACGU_num(first)
                    weight = 1.0
                    if NodeType(0, 0) in nodes_with_best_weight.get(aa, {}):
                        weight = lambda_ * math.log(edges_with_best_weight[aa][(NodeType(0, 0), first_num, n1)] / nodes_with_best_weight[aa][NodeType(0, 0)])
                    graph.add_edge(NodeType(0, 0), n1, first_num, weight)
                else:
                    n1 = NodeType(1, i-1)

                last_first = first
                second_num = get_ACGU_num(second)
                weight = 1.0
                if n1 in nodes_with_best_weight.get(aa, {}):
                    weight = lambda_ * math.log(edges_with_best_weight[aa][(n1, second_num, n2)] / nodes_with_best_weight[aa][n1])
                graph.add_edge(n1, n2, second_num, weight)

                for third in thirds:
                    three_nums = first + second + third
                    weight = 1.0
                    if n2 in nodes_with_best_weight.get(aa, {}):
                        weight = lambda_ * math.log(codon.get_weight(aa, three_nums) / nodes_with_best_weight[aa][n2])
                    else:
                        weight = lambda_ * math.log(codon.get_weight(aa, three_nums))
                    graph.add_edge(n2, NodeType(0, 0), get_ACGU_num(third), weight)

                i += 1
            aa_graphs[aa] = graph

    return aa_graphs

### (class) DFA ----------------------------------------------------------------
class DFA:
    def __init__(self):
        self.nodes = defaultdict(list)
        self.left_edges = defaultdict(list)
        self.right_edges = defaultdict(list)
        self.auxiliary_left_edges = defaultdict(lambda: defaultdict(list))
        self.auxiliary_right_edges = defaultdict(lambda: defaultdict(list))
        self.node_rightedge_weights = defaultdict(dict)

    def add_edge(self, n1, n2, nuc, weight=0.0):
        self.right_edges[n1].append(NodeNucWType(n2, nuc, weight))
        self.left_edges[n2].append(NodeNucWType(n1, nuc, weight))
        self.auxiliary_right_edges[n1][n2].append((nuc, weight))
        self.auxiliary_left_edges[n2][n1].append((nuc, weight))
        self.node_rightedge_weights[n1][nuc] = weight

    def add_node(self, n1):
        pos = n1.index
        self.nodes[pos].append(n1)

    def print(self, file=None):
        if file is None:
            import sys
            file = sys.stdout
            
        sorted_keys = sorted(self.nodes.keys())
        print("%", file=file)
        for key in sorted_keys:
            for node in self.nodes[key]:
                print(f"node:\t({node.index}, {node.num})", file=file)
                if node in self.right_edges:
                    for cnt, edge in enumerate(self.right_edges[node], 1):
                        print(f"R_{cnt}:\t{get_ACGU_char(edge.nuc)};{edge.weight};({edge.node.index}, {edge.node.num})", file=file)
                if node in self.left_edges:
                    for cnt, edge in enumerate(self.left_edges[node], 1):
                        print(f"L_{cnt}:\t{get_ACGU_char(edge.nuc)};{edge.weight};({edge.node.index}, {edge.node.num})", file=file)
                print("%", file=file)

### utils ----------------------------------------------------------------------
def convert_aa_to_triple(aa_seq):
    aa_to_triple = {
        'F': "Phe", 'L': "Leu", 'S': "Ser", 'Y': "Tyr", 
        'C': "Cys", 'W': "Trp", 'P': "Pro", 'H': "His", 
        'Q': "Gln", 'R': "Arg", 'I': "Ile", 'M': "Met",
        'T': "Thr", 'N': "Asn", 'K': "Lys", 'V': "Val", 
        'D': "Asp", 'E': "Glu", 'G': "Gly", 'A': "Ala", 
        '*': "STOP",
    }
    return ' '.join([aa_to_triple[aa] for aa in aa_seq])

def utr_to_aa(utr_seq, codon_table):
    utr_aa = ""

    for i in range(0, len(utr_seq), 3):
        codon = utr_seq[i:i+3]
        aa = codon_table[codon][0]
        utr_aa += aa

    return utr_aa

def get_ACGU_char(x):
    return 'A' if x == 1 else 'C' if x == 2 else 'G' if x == 3 else 'U' if x == 4 else 'X'

def get_ACGU_num(x):
    return 1 if x == 'A' else 2 if x == 'C' else 3 if x == 'G' else 4 if x == 'U' else 0