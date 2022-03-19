from classes.classes import Nucleotide


api_key = "bc40eac9be26ca5a6e911b42238d9a983008"
ICE = []
(acc_numbers, starts, stops) = ({},{},{})
offset = 20000
folder = "/Users/gera/Desktop/ICEs/tyrosine_recombinase/TR_evolution_manuscript/MSB_submission/revision/ICE_verification/level2/"
def parse_table(file):
    file = open(file, "r")
    line = file.readline()
    count = 0
    while line:
        line = line.rstrip()
        my_list = (line.split('\t'))
        if len(my_list)>3:
            ICE.append(my_list[2])
            acc_numbers[my_list[2]] = my_list[3]
            starts[my_list[2]] = my_list[4]
            if (int(starts[my_list[2]])-offset) < 0:
                starts[my_list[2]] = 0
            else:
                starts[my_list[2]] = int(starts[my_list[2]]) - offset
            stops[my_list[2]] = int(my_list[5]) + offset
            count += 1
        else:
             pass
        line = file.readline()
    file.close()
    #return all_accession_numbers, all_genera, genera_number

parse_table("input")
for i in ICE:
    print(i, acc_numbers[i], starts[i], stops[i])
    nucleotide = Nucleotide(acc_numbers[i])
    nucleotide.nuc_download(api_key, str(starts[i]), str(stops[i]), folder+i+".gb")