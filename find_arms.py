import os
from classes.classes import MappingTable, IdenticalProtein, Nucleotide
import argparse
import time
from Bio import SeqIO

# The script is designed to find all copies of MGEs. It requires the output of HTrate - a tool to detect horizontally
# transferred genes and uses a list of accession numbers as input. It outputs the MGEs nucleotide sequences and the
# encoded proteins. It also attempts to find target site duplications. Now looks for offtargets (do not cotain guaA
# notion in the file). To change that you have to modify
# if offtarget:  to
# if not offtarget
from collections import Counter


debug = True
start_time = time.time()
parser = argparse.ArgumentParser()
parser.add_argument("file", default="nothing", help="This is the input file")
parser.add_argument("--a", default=False, help="create alignment of nucleotide sequences")
parser.add_argument("--s", default=100, help="nucleotide sequences are expanded by this length for the alignment")
parser.add_argument("--api_key", default="none", help="This is your api_key to access NCBI")

right_end = "GAATGATTCCGCGT"
left_end = "TTTACAATAGAGTGGGA"
args = parser.parse_args()
filename = args.file
filename = os.path.abspath(filename)
directory = os.path.dirname(filename) + "/ip"
nuc_directory = os.path.dirname(filename) + "/nuc"
out_file = nuc_directory + "/All_out.fa"
protein_file = nuc_directory + "/proteins.fa"
limit_to_download = 100  # limits nuc sequences to download - should be fixed later

if not os.path.isdir(directory):
    print("You should run HTrate.py before!")
if not os.path.isdir(nuc_directory):
    print("creating %s folder" % nuc_directory)
    os.mkdir(nuc_directory)
print("your nuc_directory is", nuc_directory)
nuc_out = open(out_file, "w+")
p_out = open(protein_file, "w+")
api_key = args.api_key
print("input file is", filename)
print("your api_key is", api_key)
count, count_all, count_HT = 0, 0, 0
all_count, all_identical_array, all_identical_lens, all_genera, ht_genera = [], [], [], [], []

# here we will read the mapping table
mapping_table = MappingTable(filename, 0)
new_array = mapping_table.parse_mapping_table()
mt_length = len(new_array)
print("will process", mt_length, "accessions..")
unique = []
to_download = []
already_identical, dataframe_array, all_acc_numbers, processed_acc = [], [], [], []
# now we will download the identical protein report (IP file)
for acc_number in new_array:
    identical_protein = IdenticalProtein(acc_number, directory)
    if os.path.isfile(identical_protein.file) and os.path.getsize(identical_protein.file) > 0:
        processed_acc.append(acc_number)

print("Total to process: ", len(processed_acc))
print("Initial to parse: ", len(new_array))

# now parsing
for acc_number in new_array:
    count_all += 1
    # check if corresponding accession number does not exist in previously processed IdenticalProtein
    #    if acc_number not in already_identical:
    count += 1
    identical_protein = IdenticalProtein(acc_number, directory)
    all_count.append(count)
    unique.append(acc_number)
    # we open the IP file and remove all instances of the protein from the mapping table
    try:
        all_identical, genera, genera_number, all_nucs, copy_number, nuc_start, nuc_end, nuc_strand = identical_protein.parse_identical_protein()
    except:
        print(acc_number, "is damaged")
    counting = 0
    try:
        for v in all_nucs:
            counting += 1
            if counting > limit_to_download:
                print("reached ", limit_to_download, "sequences, breaking..")
                break
            left_border = 500
            right_border = 1000
            size = len(all_nucs)
            try:
                new_start, new_end = ("", "")
                if "+" in nuc_strand[v]:
                    start = int(nuc_start[v]) - left_border
                    end = int(nuc_end[v]) + right_border
                else:
                    start = int(nuc_start[v]) - right_border
                    end = int(nuc_end[v]) + left_border
                nucleotide = Nucleotide(v, nuc_directory, str(start), str(end), nuc_strand[v])
                if os.path.isfile(nucleotide.file) and os.path.getsize(nucleotide.file) > 0:
                    pass
                else:
                    if start > 0 and end > 0:
                        nucleotide.nuc_download(api_key, nucleotide.file)  # switch to default nuc_format and nucleotide.file for gb
                    else:
                        pass
                try:
                    tsd_min = 9
                    tsd_max = 10
                    tsd_distance = 100
                    print("looking at", nucleotide.file)
                    nucleotide.find_repeats(tsd_min, tsd_max, tsd_distance)  # here we identify repeats longer than tsd_limit
                except:
                    print("Error with find_repeats")
            except:
                pass
    except:
        pass
print("Total time: %s seconds" % (time.time() - start_time))
