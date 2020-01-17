import os
from classes.classes import MappingTable, IdenticalProtein
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

count, count_all, count_HT = 0, 0, 0
all_count = []
all_identical_array = []
all_identical_lens = []
all_genera = []
already_identical = []
genera_number = {}
dataframe_array = []
all_acc_numbers = []

# here we define input mapping table file, the only one we need for our analysis
filename = "/Users/gera/Desktop/ICEs/tyrosine_recombinase/epsilon_15/distribution_analysis/TR_distribution/SXT/SXT.mapping_table"
directory = os.path.dirname(filename)

# here we will read the mapping table
mapping_table = MappingTable(filename)
new_array = mapping_table.parse_mapping_table()
mt_length = len(new_array)

# here we define the plot
"""fig = plt.figure()
ax = fig.add_subplot(111)
Ln, = ax.plot([0, 2], color="black")
ax.set_xlim([0, mt_length])
plt.ion()
plt.show()
"""
# now we will download the identical protein report (IP file)
for acc_number in new_array:
    print(acc_number)
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
        #df = pd.DataFrame.from_dict(list(genera_number.items()), columns=['column1', 'column2'])
        print(len(genera))

        if len(genera) > 4:
            dataframe_array.append(genera_number)
            all_acc_numbers.append(acc_number)
            #df = pd.DataFrame(dataframe_array, index=all_acc_numbers)
            #df = pd.DataFrame(dataframe_array)
            count_HT += 1
        #print("====\n", df, "\n====")
        all_identical_lens.append(len(all_identical))
        all_genera.append(len(genera))
        already_identical = already_identical + all_identical
        new_array = [x for x in new_array if x not in all_identical]
        """
        if len(genera) > 4:
            Ln, = ax.plot([0, len(genera)], color="blue")
            # print("HT detected for:", acc_number, ":", genera_number)
            count_HT += 1
        Ln.set_ydata(all_genera)
        Ln.set_xdata(range(len(all_genera)))
        plt.ylabel("Number of genera")
        plt.xlabel("Index")
        plt.title("Number of genera in the protein family members")
        plt.pause(0.1)"""
        print(count_all, "out of", mt_length, " is processed")
    else:
        print(count_all, "out of", mt_length, " is already in identical")

print("total number of proteins is", count_all, "of them", count, "are unique")
print("calculated HT rate is", count_HT/count)

df = pd.DataFrame(dataframe_array, index=all_acc_numbers)
print(df)
df.plot.bar(stacked=True)
plt.show()