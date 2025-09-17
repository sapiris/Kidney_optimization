from utils import is_nan, legal_age, blood_2_blood, calc_match, parse_data, id_alleles
import numpy as np
import os


def find_best_match(top_x, transplants_file_path, year, loci_list, donor_type, f_pat_info):


    dict_don, dict_pat, pat_order_dict, don_order_dict, _, dict_id_alleles = parse_data( transplants_file_path, donor_type, year)

    unmatched_donors_count = 0
    size = len(pat_order_dict)

    #dict_id_alleles = id_alleles(dict_id_muugs)


    pat_order_dict = sorted(pat_order_dict.items(), key=lambda p: p[1])
    pat_order= [(pat_id[0], i+1) for i,pat_id in enumerate(pat_order_dict)]
    don_order_dict = sorted(don_order_dict.items(), key=lambda p: p[1])
    don_order = [(don_id[0], i + 1) for i, don_id in enumerate(don_order_dict)]
    total_hla_score = 0

    #f_out = open(f"res_top_k/top_{top_x}_{year}_{donor_type}_{('-').join(loci_list)}.csv", "w")
    #f_out.write("pat_id,don_id,order,old_order\n")
    order = 1
    list_change_inPlace = []
    while don_order:

        top_10 = pat_order[:top_x]

        max_hla_score, pat_idx, don_idx, old_order = -1, 0 , 0, 0

        #for don, don_param in dict_don.items():
        don = don_order[0]
        pat_param = dict_pat[top_10[0][0]]
        don_param = dict_don[don[0]]
        #Check if waited more than k-1 rounds
        if order - top_10[0][1]   > top_x -1 and pat_param["blood"] in blood_2_blood[don_param["blood"]] and \
                legal_age(don_param["age"], pat_param["age"]):
                    pat = top_10[0]
                    match_weight = calc_match(dict_id_alleles[don[0]], dict_id_alleles[pat[0]] , loci_list)

                    max_hla_score = match_weight
                    pat_idx = pat[0]
                    don_idx = don[0]
                    old_order = pat[1]
        else:

            for pat in top_10:
                pat_param = dict_pat[pat[0]]
                if pat_param["blood"] in blood_2_blood[don_param["blood"]]:
                    if legal_age(don_param["age"], pat_param["age"]):
                        match_weight = calc_match(dict_id_alleles[don[0]], dict_id_alleles[pat[0]], loci_list)
                        if match_weight > max_hla_score:
                            max_hla_score = match_weight
                            pat_idx = pat[0]
                            don_idx = don[0]
                            old_order = pat[1]

            if max_hla_score == -1: # if no one match by hla, take the first that match by blood and age
                for pat in pat_order:
                    pat_param = dict_pat[pat[0]]

                    if pat_param["blood"] in blood_2_blood[don_param["blood"]]:
                        if legal_age(don_param["age"], pat_param["age"]):
                            match_weight = calc_match(dict_id_alleles[don[0]], dict_id_alleles[pat[0]], loci_list)

                            max_hla_score = match_weight
                            pat_idx = pat[0]
                            don_idx = don[0]
                            old_order = pat[1]

                            break

        if max_hla_score != -1:
            total_hla_score += max_hla_score
            #f_out_ids.write(f"{pat_idx},{max_hla_score}\n")
            #f_out.write(f"{pat_idx},{don_idx},{order},{old_order}\n")
            f_info.write(f"{top_x},{('-').join(loci_list)},{pat_idx},{max_hla_score},{order},{old_order}\n")#k,loci,id,HLA score, new order, old order
            list_change_inPlace.append(order-old_order)
            order += 1
            #print(max_hla_score)


            #    print(f"{don[0]},{dict_pat[pat_idx]['blood']}, {dict_don[don[0]]['blood']}")
            pat_order.remove((pat_idx,old_order))
        else:
            unmatched_donors_count += 1
            total_hla_score+= 0#2*(len(loci_list))

        #    print(don)
        don_order.remove(don)

    # add all pat that not matched
    for pat in pat_order:
        f_info.write(f"{top_x},{('-').join(loci_list)},{pat[0]},0,-,{pat[1]}\n")
    #f_out.close()
    list_val = np.array(list_change_inPlace)

    std = list_val.std()
    change_in_place = list_val.mean()
    return round(total_hla_score/size, 2), size, std, change_in_place, unmatched_donors_count


"""f_out_score = open("hla_score_DQ_DR.csv", "w")

for i in range(1,11):
    score = find_best_match(i)
    f_out_score.write(f"{i},{score}\n")

f_out_score.close()"""

if __name__ == '__main__':
    if not os.path.exists("output"):
        os.makedirs("output")
    f_out = open("output/top_k.csv", "w")  #
    f_out.write("K,donor type, year,Avg. HLA score, Number of patients, Avg. change in order\n")

    transplants_file_path = "input_transplants.csv"



    year = 1900
    list_res_plot = {}
    list_std_plot = {}

    #check if dir exist, if no create it
    if not os.path.exists("res_top_k"):
        os.makedirs("res_top_k")

    #list_res_plot[key] = []
    #list_std_plot[key] = []
    for donor_type in [ "DDRT"]:  #choose among DDRT or LRD+DDRT
        #list_std_plot[donor_type] = []
        for k in range(0,20):
            list_res = [f"Best from {k + 1} , {donor_type} , {year}"]
            f_info = open(f"res_top_k/hla_score_per_id_{donor_type}_k={k+1}.csv", "w")
            f_info.write("k,loci,id,HLA score, new order, old order\n")
            for i, loci_list in enumerate([["A", "B", "DRB1"]]):

                #f_out = open(f"best_from_top/top10_5loci", "w")
                m = find_best_match(k+1, transplants_file_path,  year, loci_list, donor_type, f_info)
                list_res.append(str(m[0]))
                list_res.append(str(m[1]))
                list_res.append(str(m[3]))
                #list_res.append(str(m[3]))

                #print("##")
                #list_res_plot[donor_type].append(m)
                #list_std_plot[donor_type].append(std)
            list_res = (' , '). join(list_res)
            print(list_res)
            f_out.write(list_res + "\n")



