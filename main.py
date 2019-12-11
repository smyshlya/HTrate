import os
from classes.classes import MappingTable, IdenticalProtein
import matplotlib.pyplot as plt

filename = "/Users/gera/Desktop/ICEs/tyrosine_recombinase/epsilon_15/distribution_analysis/HT/Ecoli_HT.mapping_table"
directory = os.path.dirname(filename)

#here we will read the mapping table
mapping_table = MappingTable(filename)
new_array = mapping_table.parse_mapping_table()
count = 0
all_count = []
all_identical_array = []
all_identical_lens = []
all_genera = []
#now we will download the identical protein report (IP file)
for acc_number in new_array:
    count += 1
    identical_protein = IdenticalProtein(acc_number, directory)
    if os.path.isfile(identical_protein.file):
        pass
        #print("file ", identical_protein.file, " already exists\n")
    else:
        identical_protein.download()
    all_count.append(count)
#we open the IP file and remove all instances of the protein from the mapping table
    all_identical, genera = identical_protein.all_accession_numbers_and_genera()
    all_identical_lens.append(len(all_identical))
    all_genera.append(len(genera)-0)
    new_array = [x for x in new_array if x not in all_identical]
   # print("Accession numbers left to process:", new_array)

#we report in how many genera the protein is found
print(all_identical_lens)

print(all_genera)
plt.hist(all_genera, color='orange')
plt.ylabel('number of genera')
plt.show()
print("done")