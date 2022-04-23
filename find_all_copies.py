import os
from classes.classes import MappingTable, IdenticalProtein, Nucleotide
import argparse
import time
from Bio import SeqIO
debug = False
start_time = time.time()
parser = argparse.ArgumentParser()
parser.add_argument("file", default="nothing", help="This is the input file")
parser.add_argument("--a", default=False, help="create alignment of nucleotide sequences")
parser.add_argument("--s", default=100, help="nucleotide sequences are expanded by this length for the alignment")
parser.add_argument("--api_key", default="none", help="This is your api_key to access NCBI")
parser.add_argument("--ht", default=2, help="threshold number of genera for HT detection; default is 2")
args = parser.parse_args()
filename = args.file
filename = os.path.abspath(filename)
directory = os.path.dirname(filename) + "/ip"
nuc_directory = os.path.dirname(filename) + "/nuc"
out_file = nuc_directory+"/All_out.fa"
protein_file = nuc_directory+"/proteins.fa"
limit_to_download = 100  # limits nuc sequences to download - should be fixed later

if not os.path.isdir(directory):
    print("creating %s folder" % directory)
    os.mkdir(directory)
if not os.path.isdir(nuc_directory):
    print("creating %s folder" % nuc_directory)
    os.mkdir(nuc_directory)
print("your nuc_directory is", nuc_directory)
nuc_out = open(out_file, "w+")
p_out = open(protein_file, "w+")
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
    #print(all_nucs)
    try:
        #print("Total number of identical sequences to process for acc num ",acc_number,"is", len(all_nucs))

        for v in all_nucs:
            check = "CP058958"
            if check in v:
                print("processing ", v)
            counting += 1
            if counting > limit_to_download:
                print("reached ", limit_to_download, "sequences, breaking..")
                break
            left_border = 1000
            right_border = 1000
            size = len(all_nucs)
            nucleotide = Nucleotide(v, nuc_directory)
            #print("working with", nucleotide)

            try:
                new_start, new_end = ("", "")
                if "+" in nuc_strand[v]:
                    start = int(nuc_start[v]) - left_border
                    end = int(nuc_end[v]) + right_border
                else:
                    start = int(nuc_start[v]) - right_border
                    end = int(nuc_end[v]) + left_border
                if os.path.isfile(nucleotide.file) and os.path.getsize(nucleotide.file) > 0:
                    pass
                    # print(nucleotide.file, "already exists")
                else:
                    if start > 0 and end > 0 :
                        #print("trying to download", nucleotide.file)
                        nucleotide.nuc_download(api_key, str(start), str(end), nuc_strand[v], nucleotide.file)
                    else:
                        pass
                        #print(v, "is a bit out of range")
                try:
                    #print("tring to find tsds for", nucleotide.file)
                    nucleotide.find_tsd()
                except:
                    pass
                    # print("cant find tsds..")
                with open(nucleotide.file) as f:
                    genera = nucleotide.get_genera()
                    if 'horraefdaede' in f.read():
                        offtarget = False
                    else:
                        offtarget = True
                if offtarget:
                    print("offtarget")
                    position_start, position_end = "string", "string"
                    for seq_record in SeqIO.parse(nucleotide.file, "gb"):

                        position_start = seq_record.seq.find("GAATGATTCCGCGT")
                        position_end = seq_record.seq.find("TTTACAATAGAGTGGGA")
                        subseq_record = seq_record.seq[int(position_start)-400 : int(position_end)+400]
                        subRecord = seq_record[int(position_start)-400: int(position_end)+400]
                        if check in v:
                            print("for ", v, " right end is", position_end)
                            print("for ", v, " left end is", position_start)
                            print("subrecord is", subRecord)
                            print("number of features: ", len(subRecord.features))
                        for feature in subRecord.features:
                            if check in v:
                                print(feature.type)

                            if "CDS" in feature.type:
                                if check in v:
                                    print("cool")
                                    print("qualifiers are", feature.qualifiers)
                                    #print(feature)
                                #if "oxa" in feature.qualifiers['product'][0]:

                                #QLI99383
                                try:
                                    if "QLI99383" in feature.qualifiers['protein_id'][0]:
                                        print('>'+feature.qualifiers['protein_id'][0]+feature.qualifiers['product'][0],"\n", feature.qualifiers['translation'][0])
                                    #p_out.write("salut")

                                    #print("genera is ", genera)
                                    p_out.write('>'+genera+"_"+feature.qualifiers['protein_id'][0]+"_"+feature.qualifiers['product'][0]+"\n"+feature.qualifiers['translation'][0]+"\n")
                                except: print("damaged feature")
                        SeqIO.write((subRecord), nucleotide.file+"_2.gb", "genbank")

                        if subseq_record:
                            nuc_out.write(">"+seq_record.id+"\n"+str(subseq_record)+"\n")
                            #p_out.write(">" + seq_record.id + "\n" + str(subseq_record) + "\n")
                        #print(subseq_record)


            except:
                pass
                #print(v, "had some issues")
            #    print(nuc_start[v], "or", nuc_end[v], " are not integers")
        max_key = max(copy_number, key=copy_number.get)
        if copy_number[max_key] > 0:
            if debug:
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