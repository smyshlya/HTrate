# HTrate

A script that predict which bacterial proteins could have been horizontally transferred from a list of protein NCBI accession numbers (provided in the inputfile)
## Usage
```bash
usage: HTrate.py [-h] [--ht HT] [--n N] [--api_key API_KEY] file

positional arguments:
  file               This is the input file

optional arguments:
  -h, --help         show this help message and exit
  --ht HT            threshold number of genera for HT detection
  --n N              number of proteins to retrieve from the mapping table; if 0 retrieves all
  --api_key API_KEY  This is your api_key to access NCBI

```
