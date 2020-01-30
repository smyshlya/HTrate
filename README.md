# HTrate

A script that predicts, which bacterial proteins could have been horizontally transferred from a list of protein NCBI accession numbers (provided in the input file)
## Usage
```bash
usage: HTrate.py [-h] [--ht HT] [--n N] [--api_key API_KEY] file

positional arguments:
  file               input file with the list of NCBI protein accession numbers 

optional arguments:
  -h, --help         show this help message and exit
  --ht HT            threshold number of genera for HT detection
  --n N              number of proteins to retrieve from the input file; if 0 retrieves all
  --api_key API_KEY  this is your api_key to access NCBI (see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/)

```
