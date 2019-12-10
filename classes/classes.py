import os


class MappingTable:
    def __init__(self, filename):
        self.filename = filename

    def parse_mapping_table(self):
        print("uploading", self.filename,"\n")
        my_array = []
        count = 0
        file = open(self.filename, "r")
        line = file.readline()
        while line:
            line = line.rstrip()
            #print(line)
            my_array.append(line)
            count += 1
            if count > 1000:
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

    def download(self):
        api_key = "bc40eac9be26ca5a6e911b42238d9a983008"
        request = "efetch -db protein -id " + self.accession_number + " -format ipg -api_key " + api_key + " > " + self.file
        os.system(request)
        #print("requesting:", request)

    def all_accession_numbers_and_genera(self):
        file = open(self.file, "r")
        line = file.readline()
        count = 0
        all_accession_numbers = []
        all_genera = []
        while line:
            line = line.rstrip()
            my_list = (line.split('\t'))
            specie = (my_list[8])
            genera_array = (specie.split(' '))
            genera = genera_array[0]
            #print(genera)
            if genera in all_genera:
                pass
            else:
                all_genera.append(genera)
            all_accession_numbers.append(my_list[6])
            count += 1
            '''
            if count > 100:
                break
            '''
            line = file.readline()

        file.close()
        print("All genera: ", all_genera, " from ", count, " accessions")
        return all_accession_numbers
