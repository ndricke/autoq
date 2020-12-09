import os, sys

func_library = [ "C", "N", "O", "F", "Cl", "Br", # basis set doesn't include I
                "cyanide", "methylamine", # "NC",
                "dicyanamide", "nitroso",
                "OC", "carboxyl","C=O",
                "trifluoromethyl"]

Hs_dict = {
    "mepyr": [18,19,24],
    "tetry": [15,24,25,28,29,30,31,32],
    "tetrid": [30,31,34,35,36,37,38],
    }


func_file = sys.argv[1]
catalyst = "mepyr"
func_group_list = func_library  # in case we want to take a subset
func_index_list = Hs_dict[catalyst]

n = 1
func_file_list = []
for func_index_set in func_index_list:
    for func_group_set in func_group_list:

        file_name = catalyst + "_func" + str(n)
        comment = ' '.join([str(n), catalyst, str([func_index_set]), str([func_group_set]),'\n'])
        func_file_list.append(comment)
        n += 1

with open(func_file, 'w') as f:
    for item in func_file_list:
        f.write(item)
