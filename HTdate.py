import os
from classes.classes import IdenticalProtein, ProteinInstance, BioSample
import matplotlib.pyplot as plt
import pandas as pd


directory = "/Users/gera/PycharmProjects/HTrate/ip"
acc_number = "AVD07301.1"
an_folder = "/Users/gera/PycharmProjects/HTrate/an"
biosample_folder = "/Users/gera/PycharmProjects/HTrate/biosample"
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
df = pd.DataFrame()
graph = pd.DataFrame()
count=1
count2=1
for accession_number in all_accession_numbers:
    print(accession_number)
    unique_protein = ProteinInstance(accession_number, an_folder)
    count2 += 1
    if unique_protein.exists():
        print("already exists")
    else:
        unique_protein.download(api_key)
    bs_id = unique_protein.get_biosample()
    if bs_id:
        an_with_dates.append(accession_number)
        bs = BioSample(bs_id, biosample_folder)
        if bs.exists():
            print("Biosample alrady downloaded")
        else:
            bs.download(api_key)
        info = bs.get_info()
        df = df.append(info, ignore_index=True)
        count += 1
    else:
        print("Biosample does not exist for "+accession_number)
    df.info()
    if not count % 50:
        plt.close()
        df.info()
        bs.plot_info(df,'collection_date')
        """plt.figure(figsize=(20, 8))
        #df['collection_date'] = pd.to_datetime(df['collection_date'], errors='coerce')
        graph = df['location'].groupby(df['location']).count()
        graph.sort_values(axis=0, ascending=False, inplace=True)
        graph = graph.head(n=30)
        graph.plot(kind="bar")
        plt.draw()
        plt.pause(0.001)"""
    print("count"+str(count)+" and "+str(count2)+" out of "+str(len(all_accession_numbers)))
bs.plot_info(df,'collection_date')
plt.draw()
