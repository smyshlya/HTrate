import os
import subprocess
import matplotlib.pyplot as plt


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
            line = file.readline()
        file.close()
        return all_accession_numbers, all_genera, genera_number


class ProteinInstance:
    def __init__(self, accession_number, folder):
        self.an = accession_number
        self.file = folder + "/" + self.an + ".txt"

    def download(self, api_key):
        # api_key = "bc40eac9be26ca5a6e911b42238d9a983008"
        request = "efetch -db protein -id " + self.an + " -format gb -api_key " + api_key + " > " + self.file
        os.system(request)

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
                    print("found "+biosample)
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

    def __init__(self, biosample_id, folder):
        self.an = biosample_id
        self.file = folder + "/" + self.an + ".biosample"

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
        file = open(self.file, "r")
        info = {'location': 'NA', 'isolation_source': 'NA', 'collection_date': 'NA', 'sample_type': 'NA'}
        line = file.readline()
        while line:
            line = line.rstrip()
            line = line.replace("\"", "")
            if "/geographic location" in line:
                my_list = line.split('=')
                info['location'] = my_list[1]
            elif "/isolation source" in line:
                my_list = line.split('=')
                info['isolation_source'] = my_list[1]
            elif "/collection date" in line:
                my_list = line.split('=')
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

    def plot_info(self, df, value):
        plt.figure(figsize=(20, 8))
        graph = df[value].groupby(df[value]).count()
        graph.sort_values(axis=0, ascending=False, inplace=True)
        # df['collection_date'] = pd.to_datetime(df['collection_date'], errors='coerce')
        graph = graph.head(n=30)
        graph.plot(kind="bar")
        plt.draw()
        plt.pause(0.001)
