import os
from classes.classes import MappingTable, IdenticalProtein, Nucleotide
import argparse
import time
from Bio import SeqIO

start_time = time.time()
parser = argparse.ArgumentParser()
parser.add_argument("file", default="nothing", help="This is the input file")
parser.add_argument("--a", default=False, help="create alignment of nucleotide sequences")
parser.add_argument("--s", default=100, help="nucleotide sequneces are expanded by this length for the alignment")
parser.add_argument("--api_key", default="none", help="This is your api_key to access NCBI")
parser.add_argument("--ht", default=2, help="threshold number of genera for HT detection; default is 2")
args = parser.parse_args()
filename = args.file
filename = os.path.abspath(filename)
directory = os.path.dirname(filename) + "/ip"
nuc_directory = os.path.dirname(filename) + "/nuc"
out_file = nuc_directory+"/All_out.fa"

if not os.path.isdir(directory):
    print("creating %s folder" % directory)
    os.mkdir(directory)
if not os.path.isdir(nuc_directory):
    print("creating %s folder" % nuc_directory)
    os.mkdir(nuc_directory)
print("your nuc_directory is", nuc_directory)
nuc_out = open(out_file, "w+")
api_key = args.api_key
HT_threshold = int(args.ht)  # define the criteria for Horizontal Transfer
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
already_identical, dataframe_array, all_acc_numbers = [], [], []
# now we will download the identical protein report (IP file)
for acc_number in new_array:
    count_all += 1
    # check if corresponding accession number does not exist in previously processed IdenticalProtein
    if acc_number not in already_identical:
        count += 1
        identical_protein = IdenticalProtein(acc_number, directory)
        if os.path.isfile(identical_protein.file) and os.path.getsize(identical_protein.file) > 0:
            pass
        else:
            print(count)
            to_download.append(acc_number)
print("Total to download: ", len(to_download))
print("Initial to parse: ", len(new_array))

for accession_number in to_download:
    print("downloading " + accession_number)
    unique_protein = IdenticalProtein(accession_number, directory)
    unique_protein.download(api_key)

# now parsing
fasta_sequences = ""
for acc_number in new_array:
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
        all_identical, genera, genera_number, all_nucs, copy_number, nuc_start, nuc_end, nuc_strand = identical_protein.parse_identical_protein()
    except:
        print(acc_number, "is damaged")
    #now import nuc files
    counting = 0
    try:
        for v in all_nucs:
            counting += 1
            left_border = 60000
            right_border = 1000
            size = len(all_nucs)
            print(counting, "out of", size, "; nucleotide accession number is ", v, ", start is ", nuc_start[v], ", end is ", nuc_end[v], ", strand is ", nuc_strand[v])
            nucleotide = Nucleotide(v, nuc_directory)
            try:
                new_start, new_end = ("", "")
                if "+" in nuc_strand[v]:
                    print("strand is plus")
                    start = int(nuc_start[v]) - left_border
                    end = int(nuc_end[v]) + right_border
                else:
                    print("strand is minus")
                    start = int(nuc_start[v]) - right_border
                    end = int(nuc_end[v]) + left_border
                print("start and end are", start, end)
                if os.path.isfile(nucleotide.file) and os.path.getsize(nucleotide.file) > 0:
                    pass
                else:
                    if start > 0 and end > 0:
                        print("trying")
                        nucleotide.nuc_download(api_key, str(start), str(end), nuc_strand[v], nucleotide.file)
                    else:
                        print("the nucleotide sequence", v, "is too short")
                #here we check if GMP synthase is in the downloaded
                with open(nucleotide.file) as f:
                    #if 'GMP synthase' in f.read():
                    if 'horraefdaede' in f.read():
                        offtarget = False
                    else:
                        offtarget = True
                if  offtarget:
                    for seq_record in SeqIO.parse(nucleotide.file, "gb"):
                        position = seq_record.seq.find("CATTAGCGCAAGG")
                        subseq_record = seq_record.seq[int(position)-400 : int(position)+10]
                        if subseq_record:
                            nuc_out.write(">"+seq_record.id+"\n"+str(subseq_record)+"\n")
                        print(subseq_record)


            except:
                print(nuc_start[v], "or", nuc_end[v], " are not integers")
        max_key = max(copy_number, key=copy_number.get)
        if copy_number[max_key] > 0:
            print("for " + acc_number + " the maximum key is " + max_key + " : " + str(
                copy_number[max_key]) + ", genera are: " + str(genera))
    except:
        print("couldn't identify max for " + acc_number)

#here we create file with only nr sequences
seen = []
records = []
for record in SeqIO.parse(out_file, "fasta"):
    if str(record.seq) not in seen:
        seen.append(str(record.seq))
        records.append(record)


#writing to a fasta file
SeqIO.write(records, out_file+"_nr.fa", "fasta")




print("Total time: %s seconds" % (time.time() - start_time))