# https://newarkcaptain.com/how-to-retrieve-ncbi-geo-information-using-apis-part1/
# https://newarkcaptain.com/how-to-retrieve-ncbi-geo-information-using-apis-part2/

# a naive and non-comprehensive attempt on extracting sample metadata columns
# basically GSE -> all GSMs -> download and parse *_family.soft.gz per GSM

import csv
import re
from collections import Counter

import GEOparse
from Bio import Entrez

# required by NCBI
Entrez.email: str = "emre.etlioglu@gmail.com"

# gse_query: str = '((("Homo sapiens"[Organism]) AND gse[Entry Type]) AND neurodegeneration[Title]) AND ("2022/01/01"[Publication Date] : "2022/12/31"[Publication Date])'


gse_query: str = '((neurodegeneration) AND "gse"[Entry Type]) AND ("2022/01/01"[Publication Date] : "2022/12/31"[Publication Date]) '

gsm_query: str = '((gsm[Entry Type]) AND neurodegeneration[Title]) AND ("2022/12/31"[Publication Date] : "2022/03/01"[Publication Date])'


def generate_and_read_handle(query: str):
    handle = Entrez.esearch(db="gds", term=query, retmax=500, usehistory=True)
    result = Entrez.read(handle)
    handle.close()
    return result


result_gse = generate_and_read_handle(gse_query)
result_gsm = generate_and_read_handle(gsm_query)


print(f'This many publications/GSEs: {result_gse["Count"]}')
print(f'This many samples/GSMs: {result_gsm["Count"]}')

# Generate GSE identifiers
uid_regex = re.compile("[1-9]+0+([1-9]+[0-9]*)")
gse_list = ["GSE" + uid_regex.match(uid).group(1) for uid in result_gse["IdList"]]


# change this part so as to include generators
colnames = []
for gse_id in gse_list:
    gse = GEOparse.get_GEO(geo=gse_id, destdir="./")  # Get GSE object
    if not gse.get_accession() == "GSE197860": # metadata not available
        for gsm_name, gsm_data in gse.gsms.items():
            for col_name, col_value in gsm_data.metadata.items():
                colnames.append(col_name)
    del gse # does not solve the memory leak issue, need to profile the code

test = Counter(colnames)
for key, value in test.items():
    print(f"{key}: {value}")

with open("geo-sample-metadata-columns.csv", "w") as csv_file:
    writer = csv.writer(csv_file)
    for key, value in test.items():
        writer.writerow([key, value])
