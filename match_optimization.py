
import pulp
from utils import calc_match, parse_data, id_alleles, legal_age, blood_2_blood
import os
import tracemalloc
def init_ip_problem(list_don_pat, dict_don_edjs, list_pat_edjs, dict_index_pair, don_set, dict_pair_index ):
    ip_problem = pulp.LpProblem("IP_Match_Problem", pulp.LpMaximize)

    ##variables boundaries
    var = []
    for idx in dict_index_pair:
        var.append(pulp.LpVariable('x' + str(idx), lowBound=0, upBound=1, cat='Integer'))

    del dict_index_pair
    ##objection function - max: sum (x_i * w_i)
    z = 0
    for idx, pair in enumerate(list_don_pat):
        weight = pair[2]
        z += (var[idx] * weight)
    ip_problem += z

    ##
    #all pairs of don_i are no more than 1

    for don in don_set:
        don = int(don.split("_")[0])
        sum_x_of_don_i = 0
        for pat_adj in dict_don_edjs[don]:

            sum_x_of_don_i += var[dict_pair_index[(don, pat_adj)]]
        ip_problem += sum_x_of_don_i <= 1
    del don_set
    # all pairs of don_i are no more than 1

    for pat in list_pat_edjs:
        sum_x_of_pat_i = 0
        for don_adj in pat[1]:
            sum_x_of_pat_i += var[dict_pair_index[(don_adj, pat[0])]]
        ip_problem += sum_x_of_pat_i <= 1
    del list_pat_edjs
    return ip_problem, var


def create_graph_match(loci_list, transplants_file_path,  year, donor_type):


    dict_don, dict_pat, _, _, _, dict_id_alleles = parse_data( transplants_file_path, donor_type, year)

    #dict_id_alleles = id_alleles(dict_id_muugs)

    #graph = nx.Graph()
    list_don_pat = []
    dict_don_edjs = {}
    list_pat_edjs = []
    for pat, pat_param in dict_pat.items():
        list_pat = []
        pat = int(pat.split("_")[0])
        for don, don_param in dict_don.items():
            # print("pat", pat)
            # print("don", don)

            if pat_param["blood"] in blood_2_blood[don_param["blood"]]:
                if legal_age(don_param["age"], pat_param["age"]):
                    match_weight = 2* len(loci_list) - calc_match(dict_id_alleles[don], dict_id_alleles[str(pat)+ "_pat"], loci_list)
                    don = int(don.split("_")[0])
                    match_weight = match_weight if match_weight > 0 else 0.000001
                    list_don_pat.append((don, pat, match_weight))
                    if don not in dict_don_edjs:
                        dict_don_edjs[don] = []
                    dict_don_edjs[don].append(pat)
                    list_pat.append(don)
        list_pat_edjs.append([pat, list_pat])

                    #graph.add_edge(don, pat, weight=match_weight)

    # pickle.dump(graph, open('graph_don_pat.pkl', "wb"))

    dict_don =  set(dict_don.keys())
    del dict_pat
    print("Finish build graph")
    return list_don_pat, dict_don_edjs, list_pat_edjs, dict_don

if __name__ == '__main__':

    if not os.path.exists("output"):
        os.makedirs("output")
    f_out = open("output/optimal.csv", "w")  #

    transplants_file_path = "input_transplants.csv"
    year = 2010
    for donor_type in [  "dead"]:
        list_res = [f"Optimal , {donor_type} , {year} "]
        for loci_list in [["A",  "B", "DRB1"]]:

            #tracemalloc.start()

            list_don_pat, dict_don_edjs, list_pat_edjs , set_don = create_graph_match(loci_list, transplants_file_path,  year, donor_type)

            # displaying the memory
            #print("###########")
            #print(tracemalloc.get_traced_memory())


            # stopping the library
            #tracemalloc.stop()
            ##create dict of index - pair
            dict_index_pair = {}
            dict_pair_index = {}

            idx = 0
            for pair in list_don_pat:
               dict_index_pair[idx] = (pair[0], pair[1])
               dict_pair_index[(pair[0], pair[1])] = idx
               idx += 1

            sum_a = len(set_don)
            ip_problem, var = init_ip_problem(list_don_pat, dict_don_edjs, list_pat_edjs, dict_index_pair,  set_don, dict_pair_index  )

            ip_problem.solve()

            ll = ip_problem.variables()



            res_value = pulp.value(ip_problem.objective)

            #f_out = open(f"optimal/opti10_5loci", "w")

            all_w = 0
            for i, variable in enumerate(ip_problem.variables()):
                if variable.varValue > 0:

                    #print("{} = {}".format(variable.name, variable.varValue))
                    idx =int(str(variable).replace("x", ""))
                    pair = list_don_pat[idx]
                    weight = pair[2]
                    weight = weight if weight > 0.1  else 0
                    all_w += weight
                    #f_out.write(f"{pair[1]},{weight}\n")

            #print(all_w)
            #print("num pairs", sum_a)
            #print("status", ip_problem.status)
            #print("Weight match:" , round(res_value/sum_a,2) )

            list_res.append(str(round(2* len(loci_list) - res_value/sum_a,2)))
            list_res.append(str(sum_a))
        list_res = (' , ').join(list_res)
        print(list_res)
        f_out.write(list_res + "\n")


