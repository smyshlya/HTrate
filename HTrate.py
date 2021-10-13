import os
from classes.classes import MappingTable, IdenticalProtein, ProteinInstance
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import time
import re

start_time = time.time()
parser = argparse.ArgumentParser()
parser.add_argument("file", default="nothing", help="This is the input file")
parser.add_argument("--ht", default=2, help="threshold number of genera for HT detection; default is 2")
parser.add_argument("--n", default=0, help="number of proteins to retrieve from the mapping table; "
                                           "if 0 retrieves all (default)")
parser.add_argument("--api_key", default="none", help="This is your api_key to access NCBI")
args = parser.parse_args()
filename = args.file

count, count_all, count_HT = 0, 0, 0
all_count, all_identical_array, all_identical_lens, all_genera, ht_genera = [], [], [], [], []
already_identical, dataframe_array, all_acc_numbers = [], [], []
genera_number = {}
debug = True

# here we define input parameters: mapping table file, the only one we need for our analysis,
# HT threshold and protein number

filename = os.path.abspath(filename)
directory = os.path.dirname(filename) + "/ip"
if not os.path.isdir(directory):
    print("creating %s folder" % directory)
    os.mkdir(directory)
print("your directory is", directory)
HT_threshold = int(args.ht)  # define the criteria for Horizontal Transfer
threshold = int(args.n)  # how many proteins to process, if not defined process all of them
api_key = args.api_key
print("input file is", filename)
print("HT threshold is", HT_threshold)
print("number of proteins to look at is", threshold)
print("your api_key is", api_key)

# here we will read the mapping table
mapping_table = MappingTable(filename, threshold)
new_array = mapping_table.parse_mapping_table()
mt_length = len(new_array)
print("will process", mt_length, "accessions..")

unique = []
to_download = []
derefed = []
# now we will download the identical protein report (IP file)
for acc_number in new_array:
    count_all += 1
    # check if corresponding accession number does not exist in previously processed IdenticalProtein
    if acc_number not in already_identical:
        count += 1
        identical_protein = IdenticalProtein(acc_number, directory)
        if os.path.isfile(identical_protein.file) and os.path.getsize(identical_protein.file) > 0:
            derefed.append(acc_number)
#            print("++++")
        else:
            print(count)
            request = "grep -w "+acc_number+" "+directory+"/"+"DBderef.txt"
            result234 = os.popen(request).read()
            if not result234:
                to_download.append(acc_number)
                derefed.append(acc_number)
            else:
                try:
                    print("grep worked, trying now to find original accession")
                    result = re.findall(r'(.*):.*', result234)
                    print(result)
                    derefed.append(result[0])
                except ValueError:
                    print("cant't find!!!")
    else:
        derefed.append(acc_number)
print("Total to download: ", len(to_download))
print("Initial to parse: ", len(new_array))
print("Final to parse:", len(derefed))
'''for k in derefed:
    print("derefed are: ", k)'''
#ProteinInstance.download_multiple(to_download, api_key, debug, "identical", directory) # sometimes the batch download
# does not work :( then you need to activate the next for loop. It does the same but MUCH slower
for accession_number in to_download:
    print("downloading "+accession_number)
    unique_protein = IdenticalProtein(accession_number, directory)
    unique_protein.download(api_key)
#now parsing
for acc_number in derefed:
    # print(acc_number)
    count_all += 1
    # check if corresponding accession number does not exist in previously processed IdenticalProtein
#    if acc_number not in already_identical:
    count += 1
    # print(count)
    identical_protein = IdenticalProtein(acc_number, directory)
    all_count.append(count)
    unique.append(acc_number)
#        print("unique accessions: ", len(unique))
# we open the IP file and remove all instances of the protein from the mapping table
    try:
        all_identical, genera, genera_number, copy_number = identical_protein.parse_identical_protein()
    except:
        print(acc_number, "is damaged")
    try:
        max_key = max(copy_number, key=copy_number.get)
        if copy_number[max_key] > 2:
            print("for " + acc_number + " the maximum key is " + max_key + " : " + str(copy_number[max_key]))
    except:
        print("couldn't identify max for "+acc_number)

    #for key in copy_number:
    #    print(key, ' : ', copy_number[key])

    if len(genera) > HT_threshold:
        dataframe_array.append(genera_number)
        all_acc_numbers.append(acc_number)
        ht_genera.append(len(genera))
        count_HT += 1
    all_identical_lens.append(len(all_identical))
    all_genera.append(len(genera))
    already_identical = already_identical + all_identical
    new_array = [x for x in new_array if x not in all_identical]
    print(count_all, "out of", mt_length, " is processed", end='\r')
#    else:
#        print(count_all, "out of", mt_length, " is already in identical", end='\r')
    if not count % 1000:
        plt.close()

#  now we plot
        if bool(dataframe_array):
            df = pd.DataFrame(dataframe_array, index=all_acc_numbers)
            df['Total_genomes'] = df.sum(axis=1)
            df['Total_genera'] = ht_genera
            #print("xoxoxo", len(ht_genera), ht_genera)
        #    df.sort_values(['Total_genera'], axis=0, ascending=False, inplace=True)
        #    df = df.drop(['Total_genera', 'Total_genomes'], axis=1)
        #    df = df.drop('Organism', axis=1)
            df.dropna(how='any', axis=1)
            #print(df.head)
        #    df.head(n=10).plot.bar(stacked=True)
            df['Total_genera'].plot(kind="line", rot=90)
        #    plt.xticks(new_array)
                   # place text at the end of bar (subtracting 47000 from x, and 0.1 from y to make it fit within the bar)
        #    plt.annotate("xoxoxo", xy=(3, 200), color='black')

            plt.draw()
            plt.pause(0.1)
        else:
            print("No horizontally transferred proteins found")
df = pd.DataFrame(dataframe_array, index=all_acc_numbers)
df['Total_genomes'] = df.sum(axis=1)
df['Total_genera'] = ht_genera
df.dropna(how='any', axis=1)
df['Total_genera'].plot(kind="line", rot=90)
plt.savefig("HTrate_result.svg", format = 'svg')
with open(filename+".unique", 'w') as f:
    for item in unique:
        f.write("%s\n" % item)
print("total number of proteins is", count_all, "of them", count, "are unique")
print("calculated HT rate is", count_HT/count)
df.to_csv(directory + "/out.csv")
print("Total time: %s seconds" % (time.time() - start_time))
plt.show()
