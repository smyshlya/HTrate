import os
from classes.classes import IdenticalProtein, ProteinInstance, BioSample
import matplotlib.pyplot as plt
import pandas as pd
import time


start_time = time.time()
directory = "/Users/gera/PycharmProjects/HTrate/ip"
map_table_file = open("PolA_mapping_table.txt.unique", "r")
line3 = map_table_file.readline()
all_non_identical_proteins = []
count3 = 0
while line3:
    line3 = line3.rstrip()
    all_non_identical_proteins.append(line3)
    count3 += 1
    line3 = map_table_file.readline()
#    print(count3, line3)
#acc_number_array = ["VTO26435.1"]
ip_acc_number_array = ["VTO26435.1", "AVD07301.1", "WP_001120888.1", "VUX23898.1"]
#ip_acc_number_array = all_non_identical_proteins
# "VTO26435.1" - GIsul2
# "WP_001120888.1" - IS91
# "AVD07301.1" - Tn916
# "VUX23898.1" - Integron_Arginine
# "VGK33888.1" - Integron_Proline
an_folder = "/Users/gera/PycharmProjects/HTrate/an"
biosample_folder = "/Users/gera/PycharmProjects/HTrate/biosample"
df = {}
ip_count = 0


merge = False  # treat all acc_number as one df
if merge:
    count = 1  # REMOVE LATER!!! UNMUTE count = 1 further down!!
    df["merged_accessions"] = pd.DataFrame()  # REMOVE LATER!!!


# RETRIEVING  ALL ACCESSION NUMBERS
all_acc_numbers = []
all_acc_numbers_dict = {}
for ip_acc_number in ip_acc_number_array:  # now we iterate through multiple sequences
    ip_count += 1
#    print("processing ", ip_count, ip_acc_number, " out of ", len(ip_acc_number_array))
    for folder in [an_folder, biosample_folder]:
        if not os.path.isdir(folder):  # checking if 'an' folder exists
            print("creating %s folder" % folder)
            os.mkdir(folder)
    api_key = "bc40eac9be26ca5a6e911b42238d9a983008"
    identical_protein = IdenticalProtein(ip_acc_number, directory)  # creating IP object
    red_acc_numbers, genera, genera_number = identical_protein.parse_identical_protein()
    red_acc_numbers = red_acc_numbers[1:]  # all 'redundant', so to say, accession numbers in an IP record
    all_acc_numbers_dict[ip_acc_number] = red_acc_numbers
    all_acc_numbers = all_acc_numbers + red_acc_numbers
#    print("Total number of redundnat proteins in ", ip_acc_number, " is ", len(red_acc_numbers))
for keys in all_acc_numbers_dict:
    print(keys, "length of all_acc_numbers_dict is ", len(all_acc_numbers_dict[keys]))


# DOWNLOADING ALL PROTEIN FILES
all_nonRef_acc_number = []
all_uniq_prots = []
for accession_number in all_acc_numbers:  # here we iterate through all accessions of one sequence
    unique_protein = ProteinInstance(accession_number, an_folder)
    if "nr RefSeq" not in unique_protein.type:
        all_nonRef_acc_number.append(accession_number)
        if not unique_protein.exists():  # download protein gb if doesn't exist already and if it is not WP_
            print("downloading", accession_number)
            unique_protein.download(api_key)
        all_uniq_prots.append(unique_protein)
print(all_nonRef_acc_number)
print("Total number of all accession number to look at is ", len(all_acc_numbers))
print("Total number of non-reference acession numbers is", len(all_uniq_prots))
print("Trying to download multiple...")
print(all_uniq_prots)
ProteinInstance.download_multiple(all_uniq_prots, api_key)  # trying now to do batch download


# DOWNLOADING ALL CORRESPONDING BIO SAMPLE FILES
an_with_dates = []
bs_list = []
for accession_number in all_nonRef_acc_number:
    unique_protein = ProteinInstance(accession_number, an_folder)
    bs_id = unique_protein.get_biosample()  # retrieve bio sample from gb file
    if bs_id:  # check if biosample was found in gb (does not exist for reference genomes for example)
        if bs_id not in bs_list:
            an_with_dates.append(accession_number)
            bs = BioSample(bs_id, biosample_folder, accession_number)
            bs_list.append(bs)
            if not bs.exists():  # check if biosample was already downloaded. otherwise download.
                bs.download(api_key)


print("Total number of unique biosamples is ", len(bs_list))
print("Trying to download all the biosamples again..")

# now all_nonRef_acc_number is a list of accession numbers which we want to analyse

'''
    bs_id = ""
    dates = []
    an_with_dates = []
    if merge:  # define if we treat input accession numbers as one data set or create separate data sets for each
        ip_acc_number = "merged_accessions"
    else:
        df[ip_acc_number] = pd.DataFrame()  # data frame that stores all the information on all accessions for one sequence
        count = 1
    graph = pd.DataFrame()
    count2 = 1
    for accession_number in red_acc_numbers:  # here we iterate through all accessions of one sequence
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
                df[ip_acc_number] = df[ip_acc_number].append(info, ignore_index=True)  # add a row to dataframe
    #            print("MY FRAME:", df[acc_number])
                count += 1
            if not count % 1000:
                plt.close()
                df[ip_acc_number].info()
                print(df[ip_acc_number].shape)
                print(df[ip_acc_number])
                bs.plot_info(df, 'collection_date', 'all')  # last argument should be  "all" for all the years, or a specific year
    #        print("count"+str(count)+" and "+str(count2)+" out of "+str(len(all_accession_numbers)))
            outfile = ip_acc_number + ".csv"


    #bs.plot_info(df, 'collection_date', "all")
    print("writing ", outfile)
    print("count is ", count)
    df[ip_acc_number].to_csv(outfile, index=False)
    if count > 10000:
        break  # break here

# DOWNLOADING ALL BIO SAMPLES FILES FROM NCBI
# PARSING THE DOWNLOADED FILES





bs.plot_info(df, 'collection_date', 'all')
plt.savefig("result.svg", format = 'svg')
print("Total time: %s seconds" % (time.time() - start_time))


'''