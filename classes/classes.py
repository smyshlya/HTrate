import os
import subprocess
import matplotlib.pyplot as plt
import pandas as pd
import re


class MappingTable:
    def __init__(self, filename, threshold):
        self.filename = filename
        self.threshold = threshold
    def parse_mapping_table(self):
        print("uploading", self.filename,"\n")
        my_array = []
        count = 0
        file = open(self.filename, "r")
        line = file.readline()
        while line:
            line = line.rstrip()
            my_array.append(line)
            count += 1
            if count > self.threshold and self.threshold != 0:
                break
            line = file.readline()
        file.close()
        print("total of ", count, " accessions")
        return my_array


class IdenticalProtein:
    def __init__(self, accession_number, folder):
        self.accession_number = accession_number
        self.folder = folder
        self.file = folder + "/" + self.accession_number+".ip"

    def download(self, api_key):
        # api_key = "bc40eac9be26ca5a6e911b42238d9a983008"
        request = "efetch -db protein -id " + self.accession_number + " -format ipg -api_key " + api_key + " > " + self.file
        os.system(request)

    def parse_identical_protein(self):
        file = open(self.file, "r")
        line = file.readline()
        count = 0
        all_accession_numbers = []  # LIST of all accession numbers belonging to this IP
        all_genera = []  # LIST of all genera where this IP is found
        genera_number = {}  # DICTIONARY: 'genera' -> number of instances
        while line:
            line = line.rstrip()
            my_list = (line.split('\t'))
            if len(my_list)>8:
                specie = (my_list[8])
                genera_array = (specie.split(' '))
                genera = genera_array[0]
                if genera in all_genera:
                    genera_number[genera] += 1
                else:
                    genera_number[genera] = 1
                    all_genera.append(genera)
                all_accession_numbers.append(my_list[6])
                count += 1
            else:
#                print("Ipg:", self.accession_number, "is broken")
                 pass
            line = file.readline()
        file.close()
        return all_accession_numbers, all_genera, genera_number


class ProteinInstance:
    def __init__(self, accession_number, folder):
        self.an = accession_number
        self.file = folder + "/" + self.an + ".txt"
        if "WP_" in accession_number:
            self.type = "nr RefSeq"
        else:
            self.type = "usual"

    def download(self, api_key):
        # api_key = "bc40eac9be26ca5a6e911b42238d9a983008"
        request = "efetch -db protein -id " + self.an + " -format gb -api_key " + api_key + " > " + self.file
        os.system(request)

    @staticmethod
    def download_multiple(protein_instances, api_key):
        print("downloading", len(protein_instances), " proteins")

    def get_biosample(self):
        file = open(self.file, "r")
        line = file.readline()
        biosample = False
        while line:
            line = line.rstrip()
            #print(line)
            if "BioSample" in line:
                my_list = (line.split(': '))
                if len(my_list) > 1:
                    try: biosample = my_list[1]
                    except: raise Exception('Problem with '+line)
                    #print("found "+biosample)
            line = file.readline()
        file.close()
        return biosample

    def exists(self):
        if os.path.isfile(self.file) and os.path.getsize(self.file) > 0:
            return True
        else:
            return False

    def get_neighbor(self):

        return


class BioSample:

    def __init__(self, biosample_id, folder, accession_number):
        self.an = biosample_id
        self.file = folder + "/" + self.an + ".biosample"
        self.prot_an = accession_number
    def download(self, api_key):
        if os.path.isfile(self.file) and os.path.getsize(self.file) > 0:
            print("Biosample file"+self.file+"already exists")
        else:
            print("downloading "+self.an)
            request1 = "esearch -db biosample -query " + self.an + "[accn] -api_key " + api_key+" | efetch -format txt -api_key " + api_key + " > " + self.file
            process = subprocess.Popen(request1, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                                     shell=True)
            stdout, stderr = process.communicate()

    def get_info(self):
        info = {'location': 'NA', 'isolation_source': 'NA', 'collection_date': 'NA', 'sample_type': 'NA'}
        if os.path.isfile(self.file) and os.path.getsize(self.file) > 0:
            file = open(self.file, "r")
            line = file.readline()
            while line:
                line = line.rstrip()
                line = line.replace("\"", "")
                info['biosample'] = self.an
                info['accession_number'] = self.prot_an
                if "/geographic location" in line:
                    my_list = line.split('=')
                    info['location'] = my_list[1]
                elif "/isolation source" in line:
                    my_list = line.split('=')
                    info['isolation_source'] = my_list[1]
                elif "/collection date" in line:
                    my_list = line.split('=')
                    my_list2=["",""]
                    p = re.compile(r'\d+\/\d+')
                    if p.findall(my_list[1]):  # this is because times like "1900/1952" can not be transformed to Date type,
                        # so I have to cut them to just "1952"
                        print("finding them all", p.findall(my_list[1]))
                        my_list2 = my_list[1].split('/')
                        info['collection_date'] = my_list2[1]
                    # if my_list[1]
                    else:
                        info['collection_date'] = my_list[1]

                elif " /sample type" in line:
                    my_list = line.split('=')
                    info['sample_type'] = my_list[1]
                elif " /host" in line:
                    my_list = line.split('=')
                    info['host'] = my_list[1]
                line = file.readline()
        return info
        file.close()

    def exists(self):
        if os.path.isfile(self.file) and os.path.getsize(self.file) > 0:
            return True
        else:
            return False

    def plot_info(self, df, value, input_year):
        new_df = pd.DataFrame()
        formatted_date = "collection_date_formatted"
        for key in df:
            if 'collection_date' in df[key].columns:
                print("collection date column:", df[key]['collection_date'])
                df[key][formatted_date] = pd.to_datetime(df[key]['collection_date'], errors='coerce', utc=True)
            # utc=True to avoid errors on Tz-aware
            try:
                input_year
            except:
                pass
            else:
                if input_year != "all":
                    df[key] = df[key][df[key][formatted_date].dt.year == input_year]
            if value == "collection_date":

                graph = df[key][formatted_date].groupby(df[key][formatted_date].dt.year).count()
            else:
                graph = df[key][value].groupby(df[key][value]).count()
            graph = graph.to_frame()
            new_df = pd.concat([new_df, graph], axis=1)
            new_df = new_df.rename(columns={value: key, formatted_date: key})
            if value == "collection_date":
                for year in range(1932, 2020):
                    year = float('%.1f' % (year))
                    if not year in new_df.index:
                        empty_series = pd.Series(name = year)
                        new_df = new_df.append(empty_series)
        new_df = new_df.fillna(0)
        new_df = new_df.sort_index()
        if value == "collection_date":
            new_df.iloc[0:100].plot(kind="line", figsize=(10, 4), xticks=range(1932, 2020), rot=90)
            plt.yscale("log")
        else:
            new_df['Total'] = new_df.sum(axis=1)
            new_df.sort_values(by='Total', inplace=True, ascending=False)
            new_df = new_df.drop(['Total'], axis=1)
            new_df.head(n=30).plot(kind="bar", figsize=(10, 4), rot=90, fontsize='small')
        plt.draw()
        plt.pause(0.1)
