import os
from classes.classes import IdenticalProtein, ProteinInstance, BioSample
import matplotlib.pyplot as plt
import pandas as pd
import time


debug = False
start_time = time.time()
#directory = "/Users/gera/PycharmProjects/HTrate/ip"
directory = "/Users/gera/Desktop/ICEs/tyrosine_recombinase/epsilon_15/manuscript/cryoeEM_paper/Alignments/flanks/ip"
word = "carb"  # if not empty will output BioSample accession number containing the word and the exact line with the
# word (case-insensitive)

#acc_number_array = ["VTO26435.1"]
#ip_acc_number_array = ["VTO26435.1", "AVD07301.1", "WP_001120888.1", "VUX23898.1"]
ip_acc_number_array = ["WP_000954590.1"]
#ip_acc_number_array = all_non_identical_proteins
# "VTO26435.1" - GIsul2
# "WP_001120888.1" - IS91
# "AVD07301.1" - Tn916
# "VUX23898.1" - Integron_Arginine
# "VGK33888.1" - Integron_Proline
# "WP_000954590.1" - new GIsul2
# WP_001291561.1 - new Tn916
#an_folder = "/Users/gera/PycharmProjects/HTrate/an"
an_folder = "/Users/gera/Desktop/ICEs/tyrosine_recombinase/epsilon_15/manuscript/cryoeEM_paper/Alignments/flanks/an"
#biosample_folder = "/Users/gera/PycharmProjects/HTrate/biosample"
biosample_folder = "/Users/gera/Desktop/ICEs/tyrosine_recombinase/epsilon_15/manuscript/cryoeEM_paper/Alignments/flanks/biosample"
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
    # all_accession_numbers, all_genera, genera_number, all_nucs, copy_number, nuc_start, nuc_end, nuc_strand
    red_acc_numbers, genera, genera_number, all_nucs, copy_number, nuc_start, nuc_end, nuc_strand = \
        identical_protein.parse_identical_protein()
    red_acc_numbers = red_acc_numbers[1:]  # all 'redundant', so to say, accession numbers in an IP record
    all_acc_numbers_dict[ip_acc_number] = red_acc_numbers
    all_acc_numbers = all_acc_numbers + red_acc_numbers
#    print("Total number of redundant proteins in ", ip_acc_number, " is ", len(red_acc_numbers))
for keys in all_acc_numbers_dict:
    print(keys, "length of all_acc_numbers_dict is ", len(all_acc_numbers_dict[keys]))


# DOWNLOADING ALL PROTEIN FILES
all_nonRef_acc_number = []
all_uniq_prots = []
print("checking if all proteins are downloaded")
for accession_number in all_acc_numbers:  # here we iterate through all accessions of one sequence
    unique_protein = ProteinInstance(accession_number, an_folder)
    if "nr RefSeq" not in unique_protein.type:
        if debug:
            print(accession_number, "is not Ref")
        all_nonRef_acc_number.append(accession_number)
        if not unique_protein.exists():  # download protein gb if doesn't exist already and if it is not WP_
            print("downloading", accession_number)
            unique_protein.download(api_key)
        all_uniq_prots.append(accession_number)
print("done uploading\n")
print("Total number of all accession number to look at is ", len(all_acc_numbers))
print("Total number of non-reference acession numbers is", len(all_uniq_prots))
print("Trying to download multiple...")
#print(all_uniq_prots)
#ProteinInstance.download_multiple(all_uniq_prots, api_key, debug, "not", an_folder)  # trying now to do batch download


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
            if word:
                bs.find_a_word(word)
            df[ip_acc_number] = df[ip_acc_number].append(info, ignore_index=True)  # add a row to dataframe
#            print("MY FRAME:", df[acc_number])
            count += 1
        if not count % 1000:
            plt.close()
            df[ip_acc_number].info()
            print(df[ip_acc_number].shape)
            print(df[ip_acc_number])
            bs.plot_info(df, 'host', 'all')  # last argument should be  "all" for all the years,
            # or a specific year
#        print("count"+str(count)+" and "+str(count2)+" out of "+str(len(all_accession_numbers)))
        outfile = ip_acc_number + ".csv"


#bs.plot_info(df, 'collection_date', "all")
print("writing ", outfile)
print("count is ", count)
df[ip_acc_number].to_csv(outfile, index=False)


bs.plot_info(df, 'host', 'all')
plt.savefig("result.svg", format='svg')
print("Total time: %s seconds" % (time.time() - start_time))
