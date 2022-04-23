# HTrate

Python scripts to identify HTed proteins and analyse their genetic background. HTed proteins are the proteins that are present in unrelated bacteria (same protein sequence exists in different genera) and therefore are putatively horizontally transferred. 
It takes protein NCBI accession numbers (provided in the input file), creates Identical Protein (IP) record files for each of the accession numbers 
(stored in the "/ip" folder), and outputs an out.csv file with distribution of all horizontally transferred proteins.
It also creates a distribution plot for ten most widely distributed proteins. Finally it calculates the HTrate of the dataset: number of horizontally transferred protein / total number of proteins.

Running time for a dataset with ~15K protein accession numbers is approximately 3h on a regular laptop.

## Usage
```bash
usage: HTrate.py [-h] [--ht HT] [--n N] [--api_key API_KEY] file

positional arguments:
  file               input file with the list of NCBI protein accession numbers 

optional arguments:
  -h, --help         show this help message and exit
  --ht HT            threshold number of different genera for HT detection (2 is default)
  --n N              number of proteins to retrieve from the input file; if 0 retrieves all of them (0 is default)
  --api_key API_KEY  this is your api_key to access NCBI (see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) and download IP records (none is deafult, i am not sure if that works)

```
