import os
from classes.classes import IdenticalProtein, ProteinInstance, BioSample
import matplotlib.pyplot as plt
import pandas as pd
import time


start_time = time.time()
directory = "/Users/gera/PycharmProjects/HTrate/ip"
file3 = open("PolA_mapping_table.txt.unique", "r")
line3 = file3.readline()
all_genera3 = []
count3 = 0
while line3:
    line3 = line3.rstrip()
    all_genera3.append(line3)
    count3 += 1
    line3 = file3.readline()
    print(count3, line3)
#acc_number_array = ["VTO26435.1"]
acc_number_array = ["VTO26435.1", "AVD07301.1", "WP_001120888.1", "VUX23898.1"]
#acc_number_array = all_genera3
# "VTO26435.1" - GIsul2
# "WP_001120888.1" - IS91
# "AVD07301.1" - Tn916
# "VUX23898.1" - Integron_Arginine
# "VGK33888.1" - Integron_Proline
an_folder = "/Users/gera/PycharmProjects/HTrate/an"
biosample_folder = "/Users/gera/PycharmProjects/HTrate/biosample"
df = {}
count4 = 0

merge = False  # treat all acc_number as one df

# DOWNLOADING ALL FILES FROM NCBI
if merge:
    count = 1  # REMOVE LATER!!! UNMUTE count = 1 further down!!
    df["xoxoxo"] = pd.DataFrame()  # REMOVE LATER!!! UNMUTE count = 1 further down!!

for acc_number in acc_number_array:  # now we iterate through multiple sequences
    count4 += 1
    print("processing ", count4, acc_number, " out of ", len(acc_number_array))
    for folder in [an_folder, biosample_folder]:
        if not os.path.isdir(folder):
            print("creating %s folder" % folder)
            os.mkdir(folder)
    api_key = "bc40eac9be26ca5a6e911b42238d9a983008"
    identical_protein = IdenticalProtein(acc_number, directory)
    all_accession_numbers, genera, genera_number = identical_protein.parse_identical_protein()
    all_accession_numbers = all_accession_numbers[1:]
    bs_id = ""
    dates = []
    an_with_dates = []
    if merge:
        acc_number = "xoxoxo"  # REMOVE LATER!!!
    else:
        df[acc_number] = pd.DataFrame()  # data frame that stores all the information on all accessions for one sequence
        count = 1
    graph = pd.DataFrame()
    count2 = 1
    for accession_number in all_accession_numbers:  # here we iterate through all accessions of one sequence
        unique_protein = ProteinInstance(accession_number, an_folder)
        count2 += 1
        if "nr RefSeq" not in unique_protein.type:
            if not unique_protein.exists():  # download protein gb if doesn't exist already and if it is not WP_
                unique_protein.download(api_key)
            bs_id = unique_protein.get_biosample()  # retrieve biosample from gb file
            if bs_id:  # check if biosample was found in gb (does not exist for reference genomes for example)
                an_with_dates.append(accession_number)
                bs = BioSample(bs_id, biosample_folder, accession_number)
                if not bs.exists():  # check if biosample was already downloaded. otherwise download.
                    bs.download(api_key)
                info = bs.get_info()  # retrieve information from biosample
                df[acc_number] = df[acc_number].append(info, ignore_index=True)  # add a row to dataframe
    #            print("MY FRAME:", df[acc_number])
                count += 1
            if not count % 1000:
                plt.close()
                df[acc_number].info()
                print(df[acc_number].shape)
                print(df[acc_number])
                bs.plot_info(df, 'collection_date', 'all')  # last argument should be  "all" for all the years, or a specific year
    #        print("count"+str(count)+" and "+str(count2)+" out of "+str(len(all_accession_numbers)))
            outfile = acc_number + ".csv"


    #bs.plot_info(df, 'collection_date', "all")
    print("writing ", outfile)
    print("count is ", count)
    df[acc_number].to_csv(outfile, index=False)
    if count > 10000:
        break  # break here

# PARSING THE DOWNLOADED FILES

bs.plot_info(df, 'collection_date', 'all')
plt.savefig("result.svg", format = 'svg')
print("Total time: %s seconds" % (time.time() - start_time))


