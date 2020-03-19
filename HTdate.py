import os
from classes.classes import IdenticalProtein, ProteinInstance, BioSample
import matplotlib.pyplot as plt
import pandas as pd


directory = "/Users/gera/PycharmProjects/HTrate/ip"
acc_number_array = ["VTO26435.1", "AVD07301.1", "WP_001120888.1", "VUX23898.1"]
# "VTO26435.1" - GIsul2
# "WP_001120888.1" - IS91
# "AVD07301.1" - Tn916
# "VUX23898.1" - Integron_Arginine
# "VGK33888.1" - Integron_Proline
an_folder = "/Users/gera/PycharmProjects/HTrate/an"
biosample_folder = "/Users/gera/PycharmProjects/HTrate/biosample"
df = {}
for acc_number in acc_number_array:  # now we iterate through multiple sequences
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

    df[acc_number] = pd.DataFrame()  # data frame that stores all the information on all accessions for one sequence
    graph = pd.DataFrame()
    count = 1
    count2 = 1
    for accession_number in all_accession_numbers:  # here we iterate through all accessions of one sequence
        unique_protein = ProteinInstance(accession_number, an_folder)
        count2 += 1
        if not unique_protein.exists():  # download protein gb if doesn't exist already
            unique_protein.download(api_key)
        bs_id = unique_protein.get_biosample()  # retrieve biosample from gb file
        if bs_id:  # check if biosample was found in gb (does not exist for reference genomes for example)
            an_with_dates.append(accession_number)
            bs = BioSample(bs_id, biosample_folder)
            if not bs.exists():  # check if biosample was already downloaded. otherwise download.
                bs.download(api_key)
            info = bs.get_info()  # retrieve information from biosample
            df[acc_number] = df[acc_number].append(info, ignore_index=True)  # add a row to dataframe
#            print("MY FRAME:", df[acc_number])
            count += 1
        if not count % 50:
            plt.close()
            df[acc_number].info()
            print(df[acc_number].shape)
            bs.plot_info(df, 'host', acc_number)
#        print("count"+str(count)+" and "+str(count2)+" out of "+str(len(all_accession_numbers)))

plt.savefig("result.svg", format = 'svg')