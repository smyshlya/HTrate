from classes.classes import Nucleotide

# a script to download nucleotide sequences using tab file with accession number, start and stop positions eg:
# Desulfobulbus propionicus DSM 2032	Deltaproteobacteria	ICE1_Des_pro	CP002364	98557	226409
# Legionella pneumophila	Gammaproteobacteria	pLP45 	AE017354	1355533	1401059
# Legionella pneumophila	Gammaproteobacteria	ICE-βox	AE017354	2296826	2361074
# has to be debugegd!! Nucleotide class changed
input_file = "test_input_to_download.txt"  # specify your input file
api_key = ""  # specify your api key
folder = ""  # specify folder where output will be stored
ICE = []
(acc_numbers, starts, stops) = ({}, {}, {})
offset = 20000


def parse_table(file):
    file = open(file, "r")
    line = file.readline()
    count = 0
    while line:
        line = line.rstrip()
        my_list = (line.split('\t'))
        if len(my_list) > 3:
            ICE.append(my_list[2])
            acc_numbers[my_list[2]] = my_list[3]
            starts[my_list[2]] = my_list[4]
            if (int(starts[my_list[2]]) - offset) < 0:
                starts[my_list[2]] = 0
            else:
                starts[my_list[2]] = int(starts[my_list[2]]) - offset
            stops[my_list[2]] = int(my_list[5]) + offset
            count += 1
        else:
            pass
        line = file.readline()
    file.close()
    return ICE


ICE = parse_table(input_file)
for i in ICE:
    print(i, acc_numbers[i], starts[i], stops[i])
    nucleotide = Nucleotide(acc_numbers[i])
    nucleotide.nuc_download(api_key, str(starts[i]), str(stops[i]), folder + i + ".gb")