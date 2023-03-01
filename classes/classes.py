import os
import subprocess
import matplotlib.pyplot as plt
import pandas as pd
import re
import repeatfinder as rf
from Bio import Entrez
from Bio import SeqIO
import sys

# a class to work with nucleotide sequences from Genbank
class Nucleotide:
    def __init__(self, accession_number, folder, start, end, strand):
        self.an = accession_number
        self.folder = folder
        self.start = start
        self.end = end
        self.strand = strand
        self.file = folder + "/" + accession_number + "_" + self.start + "_" + self.end + "_" + self.strand + ".gb"
        self.fasta = folder + "/" + accession_number + "_" + self.start + "_" + self.end + "_" + self.strand + ".fa"

    # this function downloads the gb file that is defined by the Nucleotide object
    def nuc_download(self, api_key, output, nuc_format='gb'):
        strand = self.strand
        start = self.start
        end = self.end
        request = "efetch -db nucleotide -id " + self.an + " -format "+nuc_format+" -api_key " + api_key + " -seq_start " + (
            start) + " -seq_stop " + (end) + " -strand " + strand + "> " + output
        os.system(request)

    # this functions gets the genus name of the genome
    def get_genera(self):
        print("getting genera")
        file = open(self.file, "r")
        line = file.readline()
        generas = "no genera"
        while line:
            line = line.rstrip()
            try:
                my_list = line.split()
            except TypeError:
                my_list = ["dummy"]
            first = ""
            try:
                first = my_list[0]
            except TypeError:
                pass
            if "ORGANISM" in first:
                try:
                    generas = my_list[1]
                    break
                except TypeError:
                    print("can't find organism")
            line = file.readline()
        file.close()
        return generas

    # this function attempts to find all repeats in the Nucleotide
    def find_repeats(self, repeat_min_length, repeat_max_length, repeat_distance):
        if os.path.isfile(self.file) and os.path.getsize(self.file) > 0:
            pass
        else:
            pass
        for seq_record in SeqIO.parse(self.file, "gb"):
            my_seq = seq_record.seq
            n_max = 0
            max_substring = []
            for i in range(0, len(my_seq)):
                for k in range(repeat_min_length, repeat_max_length + 1):
                    substring = my_seq[i:i + k]
                    n = my_seq[i:i + 2*k + repeat_distance].count(substring)
                    if n == 1:
                        continue
                    else:
                        if len(substring) == k and n > n_max:
                            n_max = n
                            max_substring = [substring]
                        elif len(substring) == k and n == n_max and (substring not in max_substring):
                            max_substring.append(substring)
            print("substring", max_substring, "is repeated", n_max, "times")


    # this function checks if the Nucleotide contains a string 'gene_name' in the file.  We use it to see whether
    # the file contains a specific gene such as 'guaA'
    def find_gene_name(self, gene_name):
        print("looking for gene name:", gene_name, "in file", self.file)
        with open(self.file) as f:
            if gene_name in f.read():
                offtarget = False
            else:
                offtarget = True
        return offtarget

# this class is probably a bit of nonsense, but I will keep it for now
class MappingTable:
    def __init__(self, filename, threshold):
        self.filename = filename
        self.threshold = threshold

    def parse_mapping_table(self):
        print("uploading", self.filename, "\n")
        my_array = []
        count = 0
        file = open(self.filename, "r")
        line = file.readline()
        while line:
            line = line.rstrip()
            my_array.append(line)
            count += 1
            if count > self.threshold != 0:
                break
            line = file.readline()
        file.close()
        print("total of ", count, " accessions")
        return my_array


# a class to work with Identical Protein Group records. See https://www.ncbi.nlm.nih.gov/ipg/docs/about/
class IdenticalProtein:
    def __init__(self, accession_number, folder):
        self.accession_number = accession_number
        self.folder = folder
        self.file = folder + "/" + self.accession_number + ".ip"

    def download(self, api_key):
        request = "efetch -db protein -id " + self.accession_number + " -format ipg -api_key " + api_key + " > " + self.file
        os.system(request)

    def parse_identical_protein(self):
        file = open(self.file, "r")
        line = file.readline()
        count = 0
        all_accession_numbers = []  # LIST of all accession numbers belonging to this IP
        all_genera = []  # LIST of all genera where this IP is found
        genera_number = {}  # DICTIONARY: 'genera' -> number of instances
        all_nucs = []  # list of all nucleotide sequences where the IP is found
        copy_number = {}  # copy number per nucleotide sequence
        nuc_start, nuc_end, nuc_strand = {}, {}, {}  # start, end and strand of the nucleotides sequence,
        # corresponding to the IP
        while line:
            line = line.rstrip()
            my_list = (line.split('\t'))
            if len(my_list) > 8:
                specie = (my_list[8])
                genera_array = (specie.split(' '))
                genera = genera_array[0]
                nuc = (my_list[2])
                if nuc in all_nucs:
                    copy_number[nuc] += 1
                else:
                    copy_number[nuc] = 1
                    all_nucs.append(nuc)
                    nuc_start[nuc] = my_list[3]  # NEED TO UPDATE FOR MULTYCOPIES!!!!
                    nuc_end[nuc] = my_list[4]
                    nuc_strand[nuc] = my_list[5]
                if genera in all_genera:
                    genera_number[genera] += 1
                else:
                    genera_number[genera] = 1
                    all_genera.append(genera)
                all_accession_numbers.append(my_list[6])
                count += 1
            else:
                pass
            line = file.readline()
        file.close()
        return all_accession_numbers, all_genera, genera_number, all_nucs, copy_number, nuc_start, nuc_end, nuc_strand


# a class to work with protein sequences from Genbank
class ProteinInstance:
    def __init__(self, accession_number, folder):
        self.an = accession_number
        self.file = folder + "/" + self.an + ".txt"
        if "WP_" in accession_number:
            self.type = "nr RefSeq"
        else:
            self.type = "usual"

    # 'download' is much slower as compared to 'download_multiple' but more stable
    def download(self, api_key):
        request = "efetch -db protein -id " + self.an + " -format gb -api_key " + api_key + " > " + self.file
        os.system(request)

    # it works great and fast, but there's a problem with some protein instances that for whatever reasons could not
    # be retrieved that way - and then the whole thing crashes.
    @staticmethod
    def download_multiple(protein_instances, api_key, debug, what_you_want, folder):  # adopted mostly
        # from https://www.biostars.org/p/66921 and https://biopython.org/docs/1.74/api/Bio.Entrez.html
        #        protein_instances = ["VTO26435.1", "AVD07301.1", "VUX23898.1"]
        print("downloading", len(protein_instances), " proteins")
        n = 50
        f = open("biosample_mapping_table.txt", "w+")
        corr = open("an_corrupted.txt", "w+")
        an_to_ref = {}
        an_to_all_an = {}
        deref_file = folder + "/DBderef.txt"
        deref = open(deref_file, "a+")
        for i in range(0, len(protein_instances), n):
            new_protein_instances = protein_instances[i:i + n]
            print("looking at seqs from ", i, " to ", i + n)
            Entrez.email = "smyshlya@embl.de"
            if debug:
                print(new_protein_instances)
            request = Entrez.epost("protein", id=",".join(new_protein_instances),
                                   email="georgy.smyshlyaev@unige.ch", api_key=api_key)
            try:
                result = Entrez.read(request, validate=False)
                print("good request is", request)
            except:
                print("got problems with these accessions")
                print("bad request is", request)
            webEnv = result["WebEnv"]
            print("new WebEnv is", webEnv)
            queryKey = result["QueryKey"]
            if debug:
                print("posting successful: ", result)
                print("webEnv is", webEnv)
                print("queryKey is ", queryKey)
            if "identical" in what_you_want:
                rettype = "ipg"
                retmode = "text"
            else:
                rettype = "native"  # could be invalid
                retmode = "xml"
            records_handle = Entrez.efetch(db='protein', retmax=n,
                                           webenv=webEnv, query_key=queryKey,
                                           email="smyshlya@embl.de", api_key=api_key,
                                           retmode=retmode,
                                           rettype=rettype
                                           )
            all_lines = []
            if "identical" in what_you_want:
                line = records_handle.readline()
                line = records_handle.readline()
                while line:
                    line = line.rstrip()
                    for p in new_protein_instances:  # there's a bit of a bug now for two identical proteins being in
                        # the same genome, could get problematic
                        if p in line:
                            my_list1 = line.split('\t')
                            ip_unique_number = my_list1[0]
                            try:
                                new_out.writelines(all_lines)
                            except:
                                pass
                            all_lines = []
                            new_file = folder + "/" + p + ".ip"
                            new_out = open(new_file, "w+")
                            all_lines.append(line + "\n")
                            try:
                                deref.write("%s:%s\n" % (star, an_to_all_an[star]))
                            except:
                                pass
                            star = p
                            an_to_all_an[star] = []
                        elif "WP_" in line:
                            if star in an_to_ref:
                                pass
                            else:
                                my_list = line.split('\t')
                                an_to_ref[star] = my_list[6]
                                an_to_all_an[star].append(my_list[6])
                        else:
                            if (p == star) & (ip_unique_number in line):
                                all_lines.append(line + "\n")
                                my_list = line.split('\t')
                                an_to_all_an[star].append(my_list[6])
                    try:
                        line = records_handle.readline()
                    except:
                        pass
            else:
                records = Entrez.parse(records_handle)
                for record in records:
                    print("record has following keys:", record.keys())
                    if 'GBSeq_xrefs' in record.keys():
                        xref = record['GBSeq_xrefs']
                        print("xref is ", xref)
                        for dbs in xref:
                            if 'BioSample' in dbs['GBXref_dbname']:
                                f.write("%s\t%s\n" % (record['GBSeq_locus'], dbs['GBXref_id']))
                    else:
                        print("record has no GBSeq_xrefs")
                    corr.write("some of accesion from %s to %s are corrupted\n" % (i, i + n))
            records_handle.close
        deref.close()

    def get_biosample(self):
        file = open(self.file, "r")
        line = file.readline()
        biosample = False
        while line:
            line = line.rstrip()
            # print(line)
            if "BioSample" in line:
                my_list = (line.split(': '))
                if len(my_list) > 1:
                    try:
                        biosample = my_list[1]
                    except:
                        raise Exception('Problem with ' + line)
            line = file.readline()
        file.close()
        return biosample

    def exists(self):
        if os.path.isfile(self.file) and os.path.getsize(self.file) > 0:
            return True
        else:
            return False


class BioSample:
    def __init__(self, biosample_id, folder, accession_number):
        self.an = biosample_id
        self.file = folder + "/" + self.an + ".biosample"
        self.prot_an = accession_number

    def download(self, api_key):
        if os.path.isfile(self.file) and os.path.getsize(self.file) > 0:
            print("Biosample file" + self.file + "already exists")
        else:
            print("downloading " + self.an)
            request1 = "esearch -db biosample -query " + self.an + "[accn] -api_key " + api_key + " | efetch -format txt -api_key " + api_key + " > " + self.file
            process = subprocess.Popen(request1, stdout=subprocess.PIPE, stdin=subprocess.PIPE,
                                       shell=True)
            stdout, stderr = process.communicate()

    def get_info(self):
        info = {'location': 'NA', 'isolation_source': 'NA', 'collection_date': 'NA', 'sample_type': 'NA'}
        countries = ["USA", "Canada", "France", "China", "Australia", "United Kingdom", "Norway", "Egypt", "Germany",
                     "South Korea", "Sweden", "Brazil", "Chile", "Romania", "Japan", "India", "South Africa",
                     "Colombia",
                     "Mexico", "Belgium", "Czech Republic", "Portugal", "Nigeria", "Poland", "Thailand", "Iraq",
                     "Spain",
                     "Paraguay", "Netherlands", "Switzerland", "Myanmar"]
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
                    for i_country in countries:
                        if i_country in my_list[1]:
                            info['location'] = i_country
                elif "/isolation source" in line:
                    my_list = line.split('=')
                    info['isolation_source'] = my_list[1]
                elif "/collection date" in line:
                    my_list = line.split('=')
                    my_list2 = ["", ""]
                    p = re.compile(r'\d+\/\d+')
                    if p.findall(
                            my_list[1]):  # this is because times like "1900/1952" can not be transformed to Date type,
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
                elif " /strain" in line:
                    my_list = line.split('=')
                    info['strain'] = my_list[1]
                elif "Organism" in line:
                    my_list = line.split(':')
                    info['organism'] = my_list[1]
                line = file.readline()
        return info
        file.close()

    def find_a_word(self, word):
        if os.path.isfile(self.file) and os.path.getsize(self.file) > 0:
            file = open(self.file, "r")
            line = file.readline()
            while line:
                line = line.rstrip()
                line = line.replace("\"", "")
                if re.search(word, line, re.IGNORECASE):
                    #                if word in line:
                    print(self.an, line)
                line = file.readline()

    def exists(self):
        if os.path.isfile(self.file) and os.path.getsize(self.file) > 0:
            return True
        else:
            return False

    # now plots only collection date info..
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
                        empty_series = pd.Series(name=year)
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
            new_df.head(n=100).plot(kind="bar", figsize=(10, 4), rot=90, fontsize='small')
            plt.yscale("log")
        plt.draw()
        plt.pause(0.1)
