#the script parses arg_mge.txt file from proMGE database https://promge.embl.de/download.cgi
import re

count = {}
words = ["Int_SXT", "Int_P2", "Xer", "IntKX", "Int_CTnDOT", "Int_Des", "Int_Tn916", "Int_BPP-1"]
mge_classes = ["	IS_Tn	", "	Mobility_island	", "	Cellular	", "	Phage	", "	Phage_like	", "	Conjugative_element	", "	Integron	"]
class_count = {}
ABR_count = {}
taxa = "Proteobacteria" # leave empty for the whole dataset
filename = "/Users/gera/Downloads/arg_mge.txt"
#now we find out all ABR classes
file1 = open(filename, "r")
line = file1.readline()
ABR_classes = []
while line:
    line = line.rstrip()
    try:
        my_list = line.split("\t")
        #print(ABR_classes)
        try:
            new_ABR_classes = my_list[11].split(',')
            #print(new_ABR_classes)
            ABR_classes.extend(new_ABR_classes)
            ABR_classes = list(set(ABR_classes))
        except:
            pass
    except TypeError:
        my_list = ["dummy"]
    line = file1.readline()
print ("found following ABR classes:", ABR_classes)
#ABR_classes = ["sulfonamide_antibiotic"]
total_ABR_count = {}
for word in words:
    count[word] = 0
    class_count[word] = {}
    ABR_count[word] = {}
    for cl in mge_classes:
        class_count[word][cl] = 0
        ABR_count[word][cl] = {}
    for ABR in ABR_classes:
        ABR_count[word][ABR] = 0
        total_ABR_count[ABR] = 0

file = open(filename, "r")
line = file.readline()

while line:


    if re.search(taxa, line, re.IGNORECASE):
        line = line.rstrip()
        try:
            my_list = line.split("\t")
            # print(my_list[12])
            try:
                total_ABR_classes2 = my_list[11].split(',')
                total_ABR_classes2 = list(set(total_ABR_classes2))
                # print(ABR_classes)
                for ABR in total_ABR_classes2:
                    total_ABR_count[ABR] += 1
            except:
                pass
        except TypeError:
            my_list = ["dummy"]



        for word in words:

            if re.search(word, line, re.IGNORECASE):
                count[word]+=1

                for cl in mge_classes:
                    if re.search(cl, line, re.IGNORECASE):
                        class_count[word][cl]+=1

                try:
                    my_list = line.split("\t")
                    # print(my_list[12])
                    try:
                        ABR_classes2 = my_list[11].split(',')
                        ABR_classes2 = list(set(ABR_classes2))
                        # print(ABR_classes)
                        for ABR in ABR_classes2:
                            ABR_count[word][ABR] += 1
                    except:
                        pass
                except TypeError:
                    my_list = ["dummy"]

    line = file.readline()
for word in words:
    print(word, count[word])
    for cl in class_count[word]:
        print (cl, class_count[word][cl])
    for ABR in ABR_classes:
        print(word, ABR, "\t", ABR_count[word][ABR])

for ABR in ABR_classes:
    print("Total:", ABR, "\t", total_ABR_count[ABR])