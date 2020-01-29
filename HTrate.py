import os
from classes.classes import MappingTable, IdenticalProtein
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


count, count_all, count_HT = 0, 0, 0
all_count, all_identical_array, all_identical_lens, all_genera = [], [], [], []
already_identical, dataframe_array, all_acc_numbers = [], [], []
genera_number = {}

# here we define input parameters: mapping table file, the only one we need for our analysis, HT threshold and protein number
filename = "/Users/gera/Desktop/ICEs/tyrosine_recombinase/epsilon_15/distribution_analysis/TR_distribution/SXT/SXT.mapping_table"
directory = os.path.dirname(filename)
HT_threshold = 2  # define the criteria for Horizontal Transfer
threshold = 100  # how many proteins to process, if not defined process all of them

# here we will read the mapping table
mapping_table = MappingTable(filename, threshold)
new_array = mapping_table.parse_mapping_table()
mt_length = len(new_array)
print("will process", mt_length, "accessions..")

# now we will download the identical protein report (IP file)
for acc_number in new_array:
    #print(acc_number)
    count_all += 1
    # check if corresponding accession number does not exist in previously processed IdenticalProtein
    if acc_number not in already_identical:
        count += 1
        identical_protein = IdenticalProtein(acc_number, directory)
        if os.path.isfile(identical_protein.file) and os.path.getsize(identical_protein.file) > 0:
            pass
        else:
            identical_protein.download()
            print("downloading",acc_number, count_all, "out of", mt_length)
        all_count.append(count)

# we open the IP file and remove all instances of the protein from the mapping table
        try:
            all_identical, genera, genera_number = identical_protein.parse_identical_protein()
        except:
            print(acc_number, "is damaged")
        if len(genera) > HT_threshold:
            dataframe_array.append(genera_number)
            all_acc_numbers.append(acc_number)
            count_HT += 1
        all_identical_lens.append(len(all_identical))
        all_genera.append(len(genera))
        already_identical = already_identical + all_identical
        new_array = [x for x in new_array if x not in all_identical]
        print(count_all, "out of", mt_length, " is processed", end='\r')
    else:
        print(count_all, "out of", mt_length, " is already in identical", end='\r')

print("total number of proteins is", count_all, "of them", count, "are unique")
print("calculated HT rate is", count_HT/count)

#  now we plot
df = pd.DataFrame(dataframe_array, index=all_acc_numbers)
df['Total'] = df.sum(axis=1)
df.sort_values(['Total'], axis=0, ascending=False, inplace=True)
df = df.drop('Total', axis=1)
df.dropna(how='any', axis=1)
print(df.head)
df.plot.bar(stacked=True)
plt.show()
