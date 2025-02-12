import forgi.graph.bulge_graph as fgb
from collections import defaultdict
import networkx as nx
import matplotlib.pyplot as plt
import sys
from src.forgiworks.score_constants import *

GAP_RATIO = 25

# Define node types
node_types = {
    'f': 'start',
    't': 'end',
    's': 'stem',
    'm': 'multiloop',
    'i': 'internal_loop',
    'h': 'hairpin'
}

# return defines, merged_info
def merge_stems(graph, defines):
    parent = {}
    merged_info = defaultdict(list)  # Tracks which stems are merged into each root stem
    defines_og = defines.copy()
    graph_og = graph.copy()

    # Initialize parent for all stem nodes
    for node in graph:
        if node.startswith('s'):
            parent[node] = node
            merged_info[node].append(node)  # Each stem starts with itself
            
    def extract_index(node):
        """Extract numerical index from node name (e.g., 's8' -> 8)."""
        return int(node[1:]) if node[1:].isdigit() else float('inf')  # Handle non-numeric cases safely
    #### End of extract_index

    def find(node):
        while parent.get(node, node) != node:
            parent[node] = parent[parent[node]]  # Path compression
            node = parent[node]
        return node
    #### End of find

    def union(a, b, via):
        root_a, root_b = find(a), find(b)
        
        # Determine which root to keep based on index
        if extract_index(root_a) > extract_index(root_b):
            root_a, root_b = root_b, root_a
        
        if root_a != root_b:
            parent[root_b] = root_a
            # Merge defines positions
            defines[root_a] = sorted(defines.get(root_a, []) + defines.get(root_b, []) + defines.get(via, []))
            del defines[root_b]
            # Track merged stems
            merged_info[root_a].extend(merged_info[root_b])
            merged_info[root_a] = sorted(set(merged_info[root_a]))  # Keep unique and sorted
            del merged_info[root_b]
    #### End of union
    
    def find_endings(defines, edges):
        starting, ending = None, None
        if 'f0' in edges.keys():
            for stem in edges['f0']:
                if stem.startswith('s'):
                    starting = stem
                    break
        if 't0' in edges:
            for stem in edges['t0']:
                if stem.startswith('s'):
                    ending = stem
                    break
        if starting is None:
            if 's0' in edges.keys() and defines['s0'][0] == 1:  
                starting = 's0'
            else:
                print("No starting stem found")
                sys.exit(0)
        if ending is None:
            point = 0
            for node in edges.keys():
                if node.startswith('s') and point < defines[node][3]:
                    point = defines[node][3]
                    ending = node
        return starting, ending
    #### End of find_endings
    
    # print(find_endings(defines, graph))
    
    # Merge stems connected by interior loops
    for node, neighbors in graph.items():
        if node.startswith('i'):
            stems = [n for n in neighbors if n.startswith('s')]
            for i in range(len(stems) - 1):
                union(stems[i], stems[i + 1], node)
                
    # Special case, merge two endings if possible
    starting, ending = find_endings(defines_og, graph_og)
    for node in graph[starting]:
        if node.startswith('m'):
            if node in graph[ending]:
                union(starting, ending, node)
                break
    # End of for loop

    # Build the merged graph
    merged_graph = defaultdict(set)
    for node, neighbors in graph.items():
        if node.startswith('i'):
            continue
        root_node = find(node) if node.startswith('s') else node
        for neighbor in neighbors:
            root_neighbor = find(neighbor) if neighbor.startswith('s') else neighbor
            if root_node != root_neighbor:
                merged_graph[root_node].add(root_neighbor)
                merged_graph[root_neighbor].add(root_node)
    
    
                
    # print(merged_graph, merged_info)

    return defines, merged_info, merged_graph
#### End of merge_stems

# return filtered_data, result_graph
def organize_stems(defines, merged_info, og_graph):
    keys = defines.keys()
    for key in keys:
        if key.startswith('s'):
            temp = []
            len_list = len(defines[key])
            endval = 0
            assert len_list%2 == 0
            for x in range(len_list//2):
                if x == 0:
                    temp.extend([defines[key][x*2], defines[key][x*2+1]])
                    endval = defines[key][x*2+1]
                else:
                    if(endval + 1 == defines[key][x*2]):
                        endval = defines[key][x*2+1]
                        temp[-1] = endval
                    else:
                        temp.extend([defines[key][x*2], defines[key][x*2+1]])
                        endval = defines[key][x*2+1]
            ### end of for loop
            if len(temp) == 6:
                # ending cases
                for i in range(2):
                    temp[i] = temp[i+2]
                temp = temp[:4]
                # print(temp)
            defines[key] = temp
    filtered_data = {key: value for key, value in defines.items() if not (key.startswith('i') or key.startswith('h'))}
    
    
    keys_merged_info = merged_info.keys()
    filtered_graph = {key: value for key, value in og_graph.items() if not (key.startswith('i') or key.startswith('h'))}
    result_graph = defaultdict(set)
    
    for key, stems in filtered_graph.items():
        if key.startswith('m'):
            set_insert = set()
            for s in stems:
                if s in keys_merged_info:
                    set_insert.add(s)
                else:
                    ## stem merged to another stem
                    found = False
                    for k, v in merged_info.items():
                        for ss in v:
                            if(ss == s):
                                #print("ss: ", ss, end = " ")
                                set_insert.add(k)
                                found = True
                                break
                        if found:
                            break
                    #### End of for loop
            ## End of for loop
            # print(key, set_insert, end = " ")
            # for s in set_insert:
            #     print(s in keys_merged_info, end = " ")
            # print()
            result_graph[key] = set_insert
        elif key.startswith('s'):
            # have to erase non-key stems and i, hs
            if key not in keys_merged_info:
                continue
            else:
                filtered_set = {element for element in stems if not (element.startswith('i') or element.startswith('h') or element.startswith('s'))}
                result_graph[key] = filtered_set
        elif key.startswith('f') or key.startswith('t'):
            result_graph[key] = filtered_graph[key]
    # End of for loop
                
    return filtered_data, result_graph
#### End of organize_stems

def only_stems(defines):
    keys = defines.keys()
    tempDict = {}
    cnt = 0
    for key in keys:
        if key.startswith('s'):
            tempDict[key] = length_stem(defines[key])
            cnt+=1
    return tempDict
#### End of only_stems

def length_stem(stem):
    # [1,2,3,4]
    return min(stem[1] - stem[0], stem[3] - stem[2])
#### End of length_stem

# return path
def traverse_from_start(g, cycles):
    start_node = 'f0' if 'f0' in g else 's0'
    if start_node not in g:
        return []
    
    cycle_nodes = {node for cycle in cycles for node in cycle}
    
    path = []
    current_node = start_node
    
    while current_node not in cycle_nodes or current_node.startswith('s') or current_node == start_node:
        path.append(current_node)
        neighbors = list(g.neighbors(current_node))
        # Avoid revisiting nodes in the path
        next_nodes = [n for n in neighbors if n not in path and (n not in cycle_nodes or n.startswith('s'))]
        #print(next_nodes)
        if not next_nodes:  # If no further neighbors, stop traversal
            break
        current_node = next_nodes[0]
        
    # End of while loop
    #path.append(current_node)
    m_nodes_in_path = [node for node in path if node.startswith('m')]
    
    if len(m_nodes_in_path) > 1:
        return path
    return []
#### End of traverse_from_start

# return cycle_num, cycle_list
def get_cycles(g):
    #print(g.nodes, g.edges)
    
    cycles = list(nx.simple_cycles(nx.DiGraph(g)))
    filtered_cycles = [cycle for cycle in cycles if len(cycle) > 2]
    
    # Deduplicate cycles
    unique_cycles = {tuple(sorted(cycle)) for cycle in filtered_cycles}

    # Convert back to list of lists if needed
    unique_cycles = [list(cycle) for cycle in unique_cycles]
    
    path = traverse_from_start(g, unique_cycles)
    if path:
        # valid cycle starting from start_node
        unique_cycles.insert(0,path)
    
    return len(unique_cycles), unique_cycles
#### End of get_cycles

# return list_connecting_stems
def find_connecting_stems(edges):
    keys = edges.keys()
    list_connecting_stems = []
    for node in keys:
        if node.startswith('s'):
            if len(edges[node]) >= 3:
                # has_f = any(neigh.startswith('f') for neigh in edges[node])
                # has_t = any(neigh.startswith('t') for neigh in edges[node])
                # if not (has_f and has_t):  # Append only if it does NOT contain both 'f' and 't'
                #     list_connecting_stems.append(node)
                list_connecting_stems.append(node)
    # End of for loop
    
    return list_connecting_stems
#### End of find_connecting_stems



# return G, cycle_list, score, defines, edges, stem_lengths
def score_structure(dotbracket, name):
    
    # return max_length
    def find_max_hairpin_to_hairpin(cycle_list, stem_lengths):
        list_stem = stem_lengths.keys()
        max_length = 0
        for stem in list_stem:
            queue = []
            visited = []
            length = {key: 0 for key in list_stem}
            queue.append(stem)
            visited.append(stem)
            while queue:
                curr_stem = queue.pop(0)
                length[curr_stem] += stem_lengths[curr_stem]
                for cycle in cycle_list:
                    if curr_stem in cycle:
                        for node in cycle:
                            if node in list_stem and node not in visited:
                                queue.append(node)
                                length[node] += length[curr_stem]
                                visited.append(node)
                ## End of for loop
            ### End of while loop
            temp = max(length.values())
            # print(stem, length)
            if temp > max_length:
                max_length = temp
        return max_length 
                
    #### End of find_max_hairpin_to_hairpin
    
    # return length of longest arm
    def find_length_longest_arm(list_connecting_stems, cycle_list, defines, stem_lengths, start_cycle, target_stem):
        #print(start_cycle, target_stem, cycle_list)
        max_length = stem_lengths[target_stem]
        visited_stems = [target_stem]
        cycle_length_to_add = {}
        nowat = None
        stack = []
        num_cycles = len(cycle_list)
        if target_stem not in list_connecting_stems:
            return stem_lengths[target_stem]
        for i in range(num_cycles):
            if target_stem in cycle_list[i] and cycle_list[i] != start_cycle:
                cycle_length_to_add[i] = stem_lengths[target_stem]
                nowat = i
                break
        # End of for loop
        if not nowat:
            return stem_lengths[target_stem]
        stack.append(nowat)
        while stack:
            idx = stack.pop()
            current_cycle = cycle_list[idx]
            for node in current_cycle:
                if node.startswith('s') and node not in visited_stems and node in list_connecting_stems:
                    for i in range(num_cycles):
                        if node in cycle_list[i] and i not in cycle_length_to_add.keys():
                            cycle_length_to_add[i] = cycle_length_to_add[idx] + stem_lengths[node]
                            stack.append(i)
                            visited_stems.append(node)
                            break
                    # End of for loop
                elif node.startswith('s') and node not in visited_stems:
                    visited_stems.append(node)
                    if max_length < cycle_length_to_add[idx] + stem_lengths[node]:
                        max_length = cycle_length_to_add[idx] + stem_lengths[node]
        # End of while loop
        return max_length
    #### End of find_length_longest_arm
    
    # return score
    # calculate score of a RNA 2D structure
    def calculate_score(len_rna, cycle_list, list_connecting_stems, num_hairpin, stem_lengths, defines, edges):
        score = 0
        num_cycles = len(cycle_list)
        cycle_stem_lengths = [[] for _ in range(num_cycles)]
        cycle_stems = [[] for _ in range(num_cycles)]
        cycle_portions = [0 for _ in range(num_cycles)]
        cycle_length = [0 for _ in range(num_cycles)]
        
        #region 1: cycle parameteres
        for i in range(num_cycles):
            length = 0
            for node in cycle_list[i]:
                if node.startswith('s'):
                    if node in list_connecting_stems:
                        #print(cycle_list[i], find_length_longest_arm(list_connecting_stems, cycle_list, defines, stem_lengths, cycle_list[i], node))
                        cycle_stem_lengths[i].append(find_length_longest_arm(list_connecting_stems, cycle_list, defines, stem_lengths, cycle_list[i], node))
                        cycle_stems[i].append(node)
                        length+=find_length_longest_arm(list_connecting_stems, cycle_list, defines, stem_lengths, cycle_list[i], node)
                    elif stem_lengths[node] > len_rna * STEM_THRESHOLD_RATIO:
                        # Only long enough stems included in calculation
                        cycle_stem_lengths[i].append(stem_lengths[node])
                        cycle_stems[i].append(node)
                        length += stem_lengths[node]
                    #length += stem_lengths[node]
            ## End of for loop
            cycle_portions[i] = length / len_rna
            cycle_length[i] = sum(cycle_stem_lengths[i])
        # End of for loop
        #endregion 1
        
        # shorter max_stem_length is better
        # -10000 ~ 0, but typically -400~0
        max_h_to_h = find_max_hairpin_to_hairpin(cycle_list, stem_lengths)
        max_stem_length = max(stem_lengths.values())
        max_stem_length_ratio = max_stem_length / len_rna
        score -= (max_stem_length ** MAX_STEM_LENGTH_QUOTIENT) * MAX_STEM_LENGTH_RATE
        
        score -= (max_h_to_h ** MAX_HAIRPIN_TO_HAIRPIN_QUOTIENT) * MAX_HAIRPIN_TO_HAIRPIN_RATE
        
        # bigger num_hairpin is better, but not so critical
        # typically 10 ~ 40
        score += ((num_hairpin ** HAIRPIN_QUOTIENT) / len_rna) * HAIRPIN_RATE
        
        # multi-loop with more stems is better. multi-loop with higher portion is more critical. multi-loop with similar length of stems is better
        scores = [0 for _ in range(num_cycles)]
        
        for i in range(num_cycles):
            sum_temp = 0
            cnt_temp = 0
            for j in range(len(cycle_stem_lengths[i])):
                if cycle_stems[i][j] in list_connecting_stems:
                    continue
                sum_temp += cycle_stem_lengths[i][j]
                cnt_temp += 1
            if cnt_temp != 0:
                avg_stem_length = sum_temp / cnt_temp
            else:
                ## Connecting multiloops
                scores[i] = -1
                continue
            diff_normalized = 0 # 0~1, smaller the better
            for j in range(len(cycle_stem_lengths[i])):
                if cycle_stems[i][j] in list_connecting_stems:
                    diff_normalized += abs(avg_stem_length - cycle_stem_lengths[i][j]) / cycle_length[i] * CONNECTING_DIFF_ADV
                else:
                    diff_normalized += abs(avg_stem_length - cycle_stem_lengths[i][j]) / cycle_length[i]
            num_stem = MAX_STEM_NUM if len(cycle_stem_lengths[i]) > MAX_STEM_NUM else len(cycle_stem_lengths[i])
            scores[i] += CYCLE_SCORE_RATE * (num_stem ** CYCLE_STEM_NUM_QUOTIENT) * (cycle_portions[i] ** CYCLE_PORTION_QUOTIENT) * ((DIFF_NORM_DEFAULT - (diff_normalized+DIFF_NORM_ADD) / DIFF_NORM_DIVIDE))
        
        # multi-loop connecting stem within some range is better
        # too short stem in multi-loop should be ignored, or just increment score by small value
        
        for i in range(num_cycles):
            max_connected_score = 0
            if scores[i] == -1:
                for j in range(num_cycles):
                    for node in cycle_list[i]:
                        if node.startswith('s') and node in cycle_list[j] and i != j:
                            if max_connected_score < scores[j]:
                                max_connected_score = scores[j]
                ## End of for looop
                scores[i] = max_connected_score * CONNECTED_SCORE_RATE
            else:
                continue
        ### End of for loop
        score += sum(scores)
        
        score_min = -((len_rna * 0.5) ** MAX_HAIRPIN_TO_HAIRPIN_QUOTIENT) * MAX_HAIRPIN_TO_HAIRPIN_RATE
        score_max = 500
        
        normalized_score = (score - score_min) / (score_max - score_min) * 100
        normalized_score = 100 if normalized_score > 100 else normalized_score
        
        return normalized_score
    #### End of calculate_score
    
    len_rna = len(dotbracket)
    #print(dotbracket, flush=True)
    bg = fgb.BulgeGraph.from_dotbracket(dotbracket)
    
    defines_og = bg.defines
    edges_og = bg.edges
    
    # print("defines:")
    # print(defines)
    # print("edges:")
    # print(edges)
    
    defines, merged_info, merged_graph = merge_stems(edges_og, defines_og)
    #print(merged_graph)
    defines, edges = organize_stems(defines, merged_info, merged_graph)
    stem_lengths = only_stems(defines)
    # print("Merged and Organized defines:")
    # print(defines)
    # print("Merged Info:")
    # print(merged_info)
    # print("Merged Edges:")
    # print(edges)
    # print("s lengths:")
    # print(stem_lengths)
    
    # Initialize graph
    G = nx.Graph()

    # Add edges
    for node, neighbors in edges.items():
        for neighbor in neighbors:
            G.add_edge(node, neighbor)

    # Assign colors and sizes based on node type and length
    node_colors = []
    node_sizes = []

    for node in G.nodes():
        if node.startswith('s'):
            node_colors.append('blue')
            node_sizes.append(stem_lengths.get(node, 10) * 10)  # Scale s node sizes
        elif node.startswith('m'):
            node_colors.append('green')
            node_sizes.append(200)  # Fixed size for m nodes
        else:
            node_colors.append('red')
            node_sizes.append(200)  # Fixed size for f nodes

    # Draw the graph and save as a PDF
    # plt.figure(figsize=(12, 8))
    # pos = nx.spring_layout(G, seed=42)  # Layout for visualization
    # nx.draw(
    #     G, pos, with_labels=True, node_color=node_colors, 
    #     node_size=node_sizes, font_size=10, edge_color='gray'
    # )
    # plt.title("Graph Representation of Merged Edges")
    # plt.savefig("tmp/graph_representation_"+name+".pdf")
    
    cycle_num, cycle_list = get_cycles(G)
    
    list_connecting_stems = find_connecting_stems(edges)
    #print(len_rna, cycle_list, list_connecting_stems)
    
    score = calculate_score(len_rna, cycle_list, list_connecting_stems, len([k for k in defines_og.keys() if k.startswith('h')]), stem_lengths, defines, edges)

    # # Print the unique cycles
    # print(f"Number of unique cycles: {cycle_num}")
    # print(cycle_list)
    return G, cycle_list, score, defines, edges, stem_lengths
#### End of score_structure

# return midpoint
def get_midpoint_m(defines, edges, tofind):
    lst = defines[tofind]
    if lst:
        # content exist
        assert len(lst) > 1
        return (lst[0] + lst[1]) // 2
    # else
    s1, s2 = list(edges[tofind])
    for i in range(4):
        for j in range(4):
            if abs(defines[s1][i] - defines[s2][j]) <= 1:
                return defines[s1][i]
    return -1
#### End of get_midpoint_m

# return penalty_list
def penalties_to_test(g, cycles, defines, edges, stem_lengths, rna_length):
    penalty_list = []
    gap_size = rna_length // GAP_RATIO
    
    #print(defines, edges, stem_lengths)
    
    keys = g.nodes()
    for node in keys:
        if node.startswith('m'):
            # m cases
            #print(get_midpoint_m(defines, edges, node))
            midpoint = get_midpoint_m(defines, edges, node)
            if midpoint == -1:
                continue
            p1, p2 = midpoint - gap_size // 2, midpoint + gap_size // 2
            if not ((0 < p1 < rna_length) and (0 < p2 < rna_length)):
                continue
            else:
                penalty_list.append(str(p1)+"~"+str(p2))
                #sys.stderr.write(node + " : " + str(p1)+"~"+str(p2) + "\n")
        elif node.startswith('s'):
            #stems
            len_stem = stem_lengths[node]
            define = defines[node]
            if len_stem > gap_size * 1.1:
                margin = len_stem % gap_size if len_stem % gap_size else 2
                start1, start2 = define[0] + margin//2, define[2] + margin//2
                end1, end2 = start1 + gap_size, start2 + gap_size
                while end1 < define[1] and end1 < len_stem - 1:
                    penalty_list.append(str(start1) + "~" + str(end1))
                    #sys.stderr.write(node + " : " + str(start1)+"~"+str(end1) + "\n")
                    start1 += gap_size // 2
                    end1 = start1 + gap_size
                # End of while loop
                while end2 < define[3] and end2 < len_stem - 1:
                    penalty_list.append(str(start2) + "~" + str(end2))
                    #sys.stderr.write(node + " : " + str(start2)+"~"+str(end2) + "\n")
                    start2 += gap_size // 2
                    end2 = start2 + gap_size
                # End of while loop
    # End of for loop
    
    sorted_list = sorted(penalty_list, key=lambda x: int(x.split("~")[0]))
    return sorted_list
#### End of penalties_to_test

def main():
    dotbracket = input("dot-bracket structure: ").replace("\n", "")
    score_structure(dotbracket, "test")
##### End of main

if __name__ == '__main__':
    main()