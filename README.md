# HTrate

A Python script that idntifies proteins that are present in unrelated bacteria (same protein sequence exists in different genera). 
It takes protein NCBI accession numbers (provided in the input file), creates Identical Protein (IP) record files for each of the accession numbers 
(stored in the "/ip" folder), and outputs an out.csv file with distribution of all horizontally transferred proteins.
It also creates a plot for ten most widely distributed proteins. Finally it also calculates the HTrate.

## Usage
```bash
usage: HTrate.py [-h] [--ht HT] [--n N] [--api_key API_KEY] file

positional arguments:
  file               input file with the list of NCBI protein accession numbers 

optional arguments:
  -h, --help         show this help message and exit
  --ht HT            threshold number of different genera for HT detection
  --n N              number of proteins to retrieve from the input file; if 0 retrieves all of them
  --api_key API_KEY  this is your api_key to access NCBI (see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/) and download IP records

```
