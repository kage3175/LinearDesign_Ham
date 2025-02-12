import forgi.graph.bulge_graph as fgb
from collections import defaultdict
import subprocess

STEM_THRESHOLD = 30
GAP_PENALTY = 5

def merge_stems(graph, defines):
    parent = {}

    # Initialize parent for all stem nodes
    for node in graph:
        if node.startswith('s'):
            parent[node] = node

    def find(node):
        while parent.get(node, node) != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(a, b, via):
        root_a, root_b = find(a), find(b)
        if root_a != root_b:
            parent[root_b] = root_a
            # Merge defines positions
            defines[root_a] = sorted(defines.get(root_a, []) + defines.get(root_b, []) + defines.get(via, []))
            del defines[root_b]

    # Merge stems connected by interior loops
    # if it is hairpin stem, we should first add hairpin indexes
    for node, neighbors in graph.items():
        if node.startswith('h'):
            stems = [n for n in neighbors if n.startswith('s')]
            stem = stems[0]
            defines[stem] = sorted(defines.get(stem, []) + defines.get(node, []))
    for node, neighbors in graph.items():
        if node.startswith('i'):
            stems = [n for n in neighbors if n.startswith('s')]
            for i in range(len(stems) - 1):
                union(stems[i], stems[i + 1], node)

    # Build the merged graph
    merged_graph = defaultdict(set)
    for node, neighbors in graph.items():
        root_node = find(node) if node.startswith('s') else node
        for neighbor in neighbors:
            root_neighbor = find(neighbor) if neighbor.startswith('s') else neighbor
            if root_node != root_neighbor:
                merged_graph[root_node].add(root_neighbor)
                merged_graph[root_neighbor].add(root_node)
    
    return merged_graph, defines
#### End of merge_stems

def organize_stems(defines):
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
            if(len(temp)) == 2:
                #perfect stem
                mid = (temp[0] + temp[1]) // 2
                temp = [temp[0], mid, mid, temp[1]]
            elif (len(temp)) != 4:
                assert False
            defines[key] = temp
    return defines
#### End of organize_stems

def only_stems(defines):
    keys = defines.keys()
    tempDict = {}
    cnt = 0
    for key in keys:
        if key.startswith('s'):
            tempDict[cnt] = defines[key]
            cnt+=1
    return tempDict
#### End of only_stems

def length_stem(stem):
    # [1,2,3,4]
    return min(stem[1] - stem[0], stem[3] - stem[2])
#### End of length_stem

def parse_pstring(s):
    result = []
    for chunk in s.split(","):
        l, r = chunk.split("~")
        result.append([l,r])
    return result
#### End of parse_pstring

def how_many_split(length, starnum):
    temp = (length-GAP_PENALTY)//(GAP_PENALTY + STEM_THRESHOLD)
    return temp if temp <= starnum else starnum
#### End of how_many_split

def pntList2String(lstPnt):
    isFirst= True
    result = ""
    for p in lstPnt:
        if not isFirst:
            result += ","
        result += str(p[0]) + "~" + str(p[1])
        isFirst = False
    return result

def further(starnum):
    file_read = open("tmp/tmp_log_ld", "r")

    rnaSeq = file_read.readline().replace("\n", "")
    dotbracket = file_read.readline().replace("\n", "")
    prev_penalty_string = file_read.readline().replace("\n", "")
    assert len(prev_penalty_string) > 1 #not empty
    prev_penalty_areas = parse_pstring(prev_penalty_string)
    
    # command to run
    command = ["../LinearPartition/linearpartition", "-V", "-M"]
    
    # file
    with open("./tmp/tmp_log_dotbracket", "w") as fout:
        process = subprocess.Popen(
            command,
            stdin=subprocess.PIPE,
            stdout=fout,
            stderr=subprocess.PIPE
        )
        process.communicate(input=rnaSeq.encode())
    #### subprocess done
    
    with open("./tmp/tmp_log_dotbracket", "r") as fin:
        trash = fin.readline()
        trash = fin.readline()
        dotbracket = fin.readline().replace("\n", "").replace(" ", "")
    #### Reading dotbracket from the subprocess result done
    # for x in dotbracket:
    #     if x != '(' and x != ')' and x != '.':
    #         print(f"something wrong with x: {repr(x)}")
    #         exit(1)
    
    print("DEBUG: dotbracket by LinearPartition: " + dotbracket)

    # Parse the dot-bracket notation into a bulge graph
    bg = fgb.BulgeGraph.from_dotbracket(dotbracket)
    graph = bg.edges
    defines = bg.defines
    # print(bg.edges)
    # print(bg.defines)

    merged_graph, defines = merge_stems(graph, defines)
    # print("---- After Merging ----")
    #print("DEBUG: defines in forgi, after merging", defines, "\n")
    defines = organize_stems(defines)
    # for node, neighbors in merged_graph.items():
    #     print(f"{node}: {neighbors}")
    # print("---- After Organizing ----")
    #print("DEBUG: defines in forgi, after organizing", defines, "\n")

    defines = only_stems(defines)

    # print("---- After Stemming ----")
    #print("DEBUG: defines in forgi, after stemming", defines, "\n")
    
    # print(prev_penalty_areas)
    
    tempPntString = ""
    isFirst = True

    for elem in defines:
        if length_stem(defines[elem]) * 2 > STEM_THRESHOLD * 3 + GAP_PENALTY * 4:
            if defines[elem][1] == defines[elem][2]:
                # hairpin stem case
                # easy
                #print(defines[elem], "to split into ", end = "")
                howmany = how_many_split(length_stem(defines[elem])*2, starnum)
                size_stem = (defines[elem][3] - defines[elem][0] - GAP_PENALTY * (howmany + 1)) // howmany
                new_pnt = []
                for i in range(howmany):
                    new_pnt.append([defines[elem][0] + GAP_PENALTY*(i+1) + size_stem * i, defines[elem][0] + GAP_PENALTY*(i+1) + size_stem * (i+1)])
                #print("stem:", defines[elem], ", howmany:", howmany, "stem size:", size_stem, "new_pnt:", new_pnt, "new_pnt_string:", pntList2String(new_pnt))
                if not isFirst:
                    tempPntString += ","
                tempPntString += pntList2String(new_pnt)
                isFirst = False
                for area in prev_penalty_areas:
                    if int(area[0]) <= defines[elem][1] and  defines[elem][0] <= int(area[1]):
                        #print(area[0] + " " + area[1] + "\n", end = "")
                        break
                
    if(len(tempPntString) > 1):
        prev_penalty_string = tempPntString + "," + prev_penalty_string
        return True, prev_penalty_string
    
    return False, ""
    #print("new penalty string: " + prev_penalty_string)
        