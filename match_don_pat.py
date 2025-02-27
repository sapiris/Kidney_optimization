from utils import legal_age, blood_2_blood, calc_match, parse_data, id_alleles
import os


def pairs_sum(loci_list, transplants_file_path,  year, donor_type):

    dict_don, dict_pat, _, _, set_pairs, dict_id_alleles = parse_data(transplants_file_path,  donor_type, year)

    #dict_id_alleles = id_alleles(dict_id_muugs)
    match_all = 0
    sum_a = 0
    for pat_id, don_id in set_pairs:
        if dict_pat[pat_id]["blood"] in blood_2_blood[dict_don[don_id]["blood"]]:
            if legal_age(dict_don[don_id]["age"],  dict_pat[pat_id]["age"]):
                sum_a += 1
                weight = calc_match(dict_id_alleles[don_id], dict_id_alleles[pat_id], loci_list)
                match_all += weight
                #print(pat_id, don_id, weight)
                #f_out.write(f"{pair[0]},{weight}\n")

    #print(sum_a)
    #print("Sum all by now ", round( match_all/sum_a, 2))
    return round( match_all/sum_a, 2), sum_a


if __name__ == '__main__':
    #create path in not exist
    if not os.path.exists("output"):
        os.makedirs("output")
    f_out = open("output/current_method.csv", "w")  #
    transplants_file_path = "input_transplants.csv"

    year = 2010
    for donor_type in ["dead"]:
        list_res = [f"Current , {donor_type} , {year} "]
        for loci_list in [ ["A", "B","DRB1"]]: #
            m =   pairs_sum(loci_list, transplants_file_path, year, donor_type)
            list_res.append(str(m[0]))
            list_res.append(str(m[1]))
        list_res = (' , ').join(list_res)
        f_out.write(list_res + "\n")
        print(list_res)

