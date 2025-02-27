import pandas as pd
from utils import calc_match, parse_data, id_alleles, legal_age,  is_nan
import networkx as nx
import pickle
import pulp
from datetime import datetime
import os
def create_graph_match(loci_list, data_file_name, year, list_donor_type):
    set_pairs = set()
    data_file = pd.read_csv(data_file_name, header=0)


    """"dict_id_muugs = {}
    with open(umug_path) as muug_file:
        for line in muug_file:
            line = line.strip().split(',')
            if not line[0] in dict_id_muugs:
                dict_id_muugs[line[0]] = {}
            dict_id_muugs[line[0]][line[1]] = float(line[2])


    dict_id_alleles = {}
    for id, dict_muugs in dict_id_muugs.items():
        dict_alleles = {"A": {}, "B": {}, "C":{}, "DQB1": {},  "DRB1": {} }#"DPB1": {},
        sum_all = sum(list(dict_muugs.values()))
        for muug, freq in dict_muugs.items():
            muug = muug.split('^')
            for allele in muug:
                locus = allele.split("*")[0]
                dict_alleles[locus][allele] = dict_alleles[locus].get(allele, 0) + freq/sum_all

        id = id.replace(".0", "")
        dict_id_alleles[str(id)] = dict_alleles"""

    dict_don = {}
    dict_pat = {}
    dict_id_alleles = {}

    blood_2_blood = {"O": ["A", "B", "AB", "O"],
                     "A": ["A", "AB"],
                     "B": ["B", "AB"],
                     "AB": ["AB"]}


    pairs = 0

    dict_date_group = {}
    group = 0
    for year in range(year, 2023):#[2020,2021,2022]:
        for month in [1,5,9]:
            group += 1
            for i in range(4):
                dict_date_group[(month+i, year)] = group

    graph = nx.DiGraph()
    sum_weight_match = 0
    for row in data_file.iterrows():
        row = row[1]
        id_line = row["SN"]
        if row["Donor type"] in list_donor_type:
            if row["Sensitive"] == "No" and  row["Cross"] == "Negative":# 0 == row["Previous transplants"]: #is_nan(row["Previous transplants"]) or
                if not is_nan(row["blood type donor"]):
                    if row["blood type recipient"] in blood_2_blood[row["blood type donor"]]:
                        #if f"{id_line}_pat" in dict_id_muugs and f"{id_line}_don" in dict_id_muugs:
                            date_format = '%Y-%m-%d'
                            date_obj = datetime.strptime(row["Date"], date_format)
                            if (date_obj.month, date_obj.year ) in dict_date_group:# and dict_date_group[(row["Date"].month, row["Date"].year )] > 1 :

                                dict_alleles_pat, dict_alleles_don = {"A": {}, "B": {}, "DRB1": {}}, {"A": {}, "B": {},
                                                                                                      "DRB1": {}}
                                for locus in ["A", "B", "DRB1"]:
                                    dict_alleles_pat[locus][row[f"HLA-{locus} recipient"]] = 1
                                    dict_alleles_don[locus][row[f"HLA-{locus} donor"]] = 1
                                dict_id_alleles[f"{id_line}_pat"] = dict_alleles_pat
                                dict_id_alleles[f"{id_line}_don"] = dict_alleles_don
                                match_weight = calc_match(dict_alleles_don, dict_alleles_pat, loci_list)
                                sum_weight_match += (6 - match_weight)
                                #if match_weight == 0:
                                #   match_weight = 0.01
                                graph.add_node(id_line, pat_blood=row["blood type recipient"], don_blood=row["blood type donor"],
                                               pat_age=row["age recipient"], don_age=row["age donor"], group_date = dict_date_group[(date_obj.month, date_obj.year )])
                                #match_weight = match_weight if match_weight > 0 else 0.000001
                                graph.add_edge(id_line, id_line, weight=match_weight)
                                """set_pairs.add((f"{id_line}_pat", f"{id_line}_don"))
                                dict_pat[str(f"{id_line}_pat")] = {"blood": row["blood type"], "age": row["age"]}
                                dict_don[str(f"{id_line}_don")] = {"blood": row["blood type.1"], "age": row["age.1"]}"""
                                pairs += 1


    print("sum match", sum_weight_match/pairs)
    print("pairs", pairs)

    pat_bloods = nx.get_node_attributes(graph, "pat_blood")
    don_bloods = nx.get_node_attributes(graph, "don_blood")
    pat_ages = nx.get_node_attributes(graph, "pat_age")
    don_ages = nx.get_node_attributes(graph, "don_age")
    dates_group = nx.get_node_attributes(graph, "group_date")

    for pat in graph.nodes():
        for don in graph.nodes():
            if pat != don:
                if pat_bloods[pat] in blood_2_blood[don_bloods[don]]:
                    if legal_age(don_ages[don], pat_ages[pat]):
                        if dates_group[don] == dates_group[pat]:
                            match_weight = calc_match(dict_id_alleles[str(f"{don}_don")],
                                              dict_id_alleles[str(f"{pat}_pat")], loci_list)
                            if match_weight > graph.get_edge_data(pat, pat)["weight"]:
                                graph.add_edge(don, pat, weight=match_weight)


    print("Finish build graph")
    return graph, sum_weight_match/pairs, pairs#, set(dict_don.keys()), set(dict_pat.keys()), set_pairs

def find_cycles(graph, k, cycle_type, dict_cycles, idx, dict_node_cycles ):
    #cycles_list = pickle.load(open('all_cycles_related.pkl', "rb"))
    cycles_list = list(nx.simple_cycles(graph))

    for cycle in cycles_list:
        if len(cycle) <= k:
            cycle_len = len(cycle)
            match_score = 0
            for i in range(cycle_len):
                match_score += graph.get_edge_data(cycle[i], cycle[(i+1)%cycle_len])["weight"]
                if cycle[i] not in dict_node_cycles: dict_node_cycles[cycle[i]] = []
                dict_node_cycles[cycle[i]].append(idx)
            dict_cycles[idx] = [cycle, match_score]
            idx += 1
    """pickle.dump(cycles_list, open(f'pkl/all_cycles_related_{cycle_type}.pkl', "wb"))
    pickle.dump(graph, open(f'pkl/graph_related_{cycle_type}.pkl', "wb"))
    pickle.dump(dict_cycles, open(f'pkl/dict_cycles_related_{cycle_type}.pkl', "wb"))
    pickle.dump(dict_node_cycles, open(f'pkl/dict_node_cycles_related_{cycle_type}.pkl', "wb"))"""

    return dict_cycles, dict_node_cycles


def init_ip_problem(graph, dict_cycles, dict_node_cycles):
    ip_problem = pulp.LpProblem("IP_Match_Problem", pulp.LpMaximize)

    ##variables boundaries and
    ##objection function init - max: sum (x_i * w_i)
    var = []
    z = 0
    for idx, info in dict_cycles.items():
        var.append(pulp.LpVariable('x' + str(idx), lowBound=0, upBound=1, cat='Integer'))
        z += (var[idx] * info[1])
    ip_problem += z

    #all cycles with pair i are no more than 1

    for node, cycles in dict_node_cycles.items():
        sum_x_of_cycle_i = 0
        for node2 in cycles:
            sum_x_of_cycle_i += var[node2]
        ip_problem += sum_x_of_cycle_i <= 1

    return ip_problem, var

#main
if __name__ == '__main__':
    if not os.path.exists("output"):
        os.makedirs("output")
    f_out = open("output/related.csv", "w")  #
    transplants_file_path = "input_transplants.csv"
    year = 2010
    for dnrtypy in [["LRD"], ["LURD"], ["LRD", "LURD"]]:
        for loci_list in [ ["A", "B","DRB1" ]]: #
            for k in [2,3,4,5,8]:
                graph_matching, sum_weights, num_pairs = create_graph_match(loci_list, transplants_file_path,  year, dnrtypy)
                dict_cycles = {}
                idx = 0
                dict_node_cycles = {}

                dict_cycles, dict_node_cycles =find_cycles(graph_matching, k, "simple" ,  dict_cycles, idx, dict_node_cycles)


                ip_problem, var = init_ip_problem(graph_matching, dict_cycles, dict_node_cycles)

                ip_problem.solve()

                ll = ip_problem.variables()

                res_value = pulp.value(ip_problem.objective)



                sum_a = 0
                all_w_cycle = 0
                all_w_orig = 0
                dict_size = {}

                for i, variable in enumerate(ip_problem.variables()):
                    if variable.varValue > 0:
                        #print("{} = {}".format(variable.name, variable.varValue))
                        cycle_info = dict_cycles[int(str(variable).replace("x", ""))]
                        size = len(cycle_info[0])
                        dict_size[size] = dict_size.get(size, 0) + 1

                        if size > 1:
                            all_w_cycle += cycle_info[1]
                            sum_a += size
                            for g in cycle_info[0]:
                                all_w_orig += graph_matching.get_edge_data(g,g)["weight"]

                #print(dict_size)
                #print(res_value/len(graph_matching.nodes()))
                #print(f"{all_w_cycle},{all_w_cycle/sum_a},{sum_a}")
                #print(f"{all_w_orig},{all_w_orig/sum_a},{sum_a}")
                #sort dict and convert to string
                sum_all =0
                dict_size = dict(sorted(dict_size.items(), key=lambda x: x[0]))
                for cycle_size, v in dict_size.items():
                    sum_all = sum_all + v*cycle_size
                if num_pairs > sum_all:
                    dict_size[1] += (num_pairs - sum_all)
                dict_size = str(dict_size).replace("{", "").replace("}", "")
                #dict_size = str(dict_size).replace("[", "").replace("]", "").replace("(", "").replace(")", "").replace(" ", "")
                #print(f"{('_').join(loci_list)} & {k} & {round( sum_weights, 2)} &{round(res_value/len(graph_matching.nodes()), 2)} & {round(all_w_orig/sum_a,2)} & {round( all_w_cycle/sum_a, 2)} & {dict_size} & {sum_all}\\\\ ")
                num_alleles = 2 * len(loci_list)
                if sum_a == 0:
                    f_out.write(f"{'_'.join(dnrtypy)}  & {('_').join(loci_list)} & {k} & {num_alleles - round(sum_weights, 2)} &{num_alleles - (round(res_value / len(graph_matching.nodes()), 2))} & -- & -- & {dict_size}  & {num_pairs}\\\\ \n")
                else:
                    f_out.write(f"{'_'.join(dnrtypy)} & {('_').join(loci_list)} & {k} & { round(num_alleles - sum_weights, 2)} &{ (round(num_alleles - (res_value / len(graph_matching.nodes())), 2))} & "
                                f"{(round(num_alleles - (all_w_orig / sum_a), 2))} & {(round(num_alleles - (all_w_cycle / sum_a), 2))} & {dict_size}  & {num_pairs}\\\\ \n")

