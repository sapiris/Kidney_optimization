import pandas as pd
from datetime import datetime
def is_nan(x):
    if x != x:
        return True
    return False


def legal_age(don_age, pat_age):

    if don_age <= 18:
        if pat_age <= 18:
            return True
        return False

    if don_age >= 65:
        if pat_age >= 60:
            return True
        return False

    return True


blood_2_blood = {"O": [ "O"],
                 "A": ["A"],
                 "B": [ "B"],
                 "AB": ["AB"]}


def calc_match(don_alleles, pat_alleles, loci_list):
    match_sum = 0

    for locus in loci_list:# ["A", "C","B","DRB1", "DQB1" ]: #don_alleles: #, "DPB1""A", "C",,"DQB1"
        for alleles_pair_don, freq_pair_don in don_alleles[locus].items():
            alleles_don = alleles_pair_don.split("+")
            for alleles_pair_pat, freq_pair_pat in pat_alleles[locus].items():
                alleles_pat = alleles_pair_pat.split("+")
                match = 0
                for allele in alleles_don:
                    if allele not in alleles_pat:#if allele in alleles_pat:
                        match += 1
                if match > 0:
                    match_sum = match_sum + match * freq_pair_don * freq_pair_pat

    return match_sum


def parse_data(transplantation_file_path,  donor_type, year = 1900):
    """dict_id_muugs = {}
    with open(hla_file_path) as muug_file: #"../belinson/output/don.belinson.umug"
        for line in muug_file:
            line = line.strip().split(',')
            if not line[0] in dict_id_muugs:
                dict_id_muugs[line[0]] = {}
            dict_id_muugs[line[0]][line[1]] = float(line[2])"""

    dict_don, dict_pat, pat_order_dict, don_order_dict  = {}, {}, {}, {}
    set_pairs = set()
    dict_id_alleles = {}
    data_file = pd.read_csv( transplantation_file_path, header=0) #"../belinson/HLA_2.xlsx"
    for row in data_file.iterrows():
        row = row[1]
        id_line = row["SN"]
        if (donor_type == "altro" and row["Donor type"] == "LURD") or (
                donor_type == "dead" and row["Donor type"] == "DDRT") or (donor_type == "altro+dead" and
                                                                                  row["Donor type"] in ["LURD", "DDRT"]):
            if row["Cross"] == "Negative":  # remove crossover transplants
                if row["Sensitive"] == "No":# 0 == row["Previous transplants"]:  # Leave only the first transplant to remove patients with antibodies. is_nan(row["Previous transplants"]) or
                    if not is_nan(row["blood type donor"]) and not is_nan(row["blood type recipient"]):
                        if row["blood type recipient"] in blood_2_blood[row["blood type donor"]]:
                            if legal_age(row["age donor"],  row["age recipient"]):
                            #if f"{id_line}_pat" in dict_id_muugs and f"{id_line}_don" in dict_id_muugs: #leave only patients with imputed HLA
                                date_format = '%Y-%m-%d'
                                date_obj = datetime.strptime(row["Date"], date_format)
                                if date_obj.year >= year:
                                    set_pairs.add((f"{id_line}_pat", f"{id_line}_don"))
                                    pat_order_dict[str(f"{id_line}_pat")] = date_obj
                                    don_order_dict[str(f"{id_line}_don")] = date_obj
                                    dict_pat[str(f"{id_line}_pat")] = {"blood": row["blood type recipient"], "age": row["age recipient"]}
                                    dict_don[str(f"{id_line}_don")] = {"blood": row["blood type donor"],
                                                                       "age": row["age donor"]}
                                    dict_alleles_pat, dict_alleles_don = {"A": {}, "B": {}, "DRB1": {}}, {"A": {}, "B": {}, "DRB1": {}}
                                    for locus in ["A", "B", "DRB1"]:
                                        dict_alleles_pat[locus][row[f"HLA-{locus} recipient"]] = 1
                                        dict_alleles_don[locus][row[f"HLA-{locus} donor"]] = 1
                                    dict_id_alleles[f"{id_line}_pat"] = dict_alleles_pat
                                    dict_id_alleles[f"{id_line}_don"] = dict_alleles_don

    return dict_don, dict_pat, pat_order_dict, don_order_dict, set_pairs, dict_id_alleles



def id_alleles(dict_id_muugs):
    dict_id_alleles = {}
    for id, dict_muugs in dict_id_muugs.items():
        dict_alleles = {"A": {}, "B": {}, "C": {}, "DQB1": {}, "DRB1": {}}  # "DPB1": {},
        sum_all = sum(list(dict_muugs.values()))
        for muug, freq in dict_muugs.items():
            muug = muug.split('^')
            for allele in muug:
                locus = allele.split("*")[0]
                if locus in dict_alleles:
                    dict_alleles[locus][allele] = dict_alleles[locus].get(allele, 0) + freq / sum_all
        id = id.replace(".0", "")
        dict_id_alleles[str(id)] = dict_alleles

    return dict_id_alleles
