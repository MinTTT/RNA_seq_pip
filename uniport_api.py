# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""
import json
# Built-in/Generic Imports
import os
import sys
# […]
#%%
# Libs
import pandas as pd
import numpy as np  # Or any other
# […]

# Own modules
import urllib.parse
import urllib.request
import json
import urllib.request
"""
the columns titles for programmatic access: https://www.uniprot.org/help/uniprotkb_column_names
Batch retrieval of entries: https://www.uniprot.org/uploadlists/
database identifier https://www.uniprot.org/help/api_idmapping

"""


# 1. Use the correct search endpoint
base_url = 'https://rest.uniprot.org/uniprotkb/search'


def query_uniprot_id(query_id):

    # query_id = "P03709"  # example query, can be a UniProt ID, gene name, etc.
    # Define search parameters
    params = {
        "query": query_id,
        "format": "json",
        "fields": "accession,protein_name,gene_names,organism_name,mass,sequence",
        "size": "500"  # recommended for performance
    }

    # 2. Encode parameters into the URL for a GET request
    query_string = urllib.parse.urlencode(params)
    full_url = f"{base_url}?{query_string}"

    # 3. Add a User-Agent (recommended by UniProt)
    req = urllib.request.Request(full_url, headers={'User-Agent': 'Python Script'})

    try:
        with urllib.request.urlopen(req) as response:
            result = response.read().decode('utf-8')
            # only get one result, if there are multiple results, you can use the "results" field in the json response to get all results
            """
            {"results":[{"entryType":"UniProtKB reviewed (Swiss-Prot)","primaryAccession":"P03709","organism":{"scientificName":"Escherichia phage lambda","commonName":"Bacteriophage lambda","taxonId":2681611,"lineage":["Viruses","Duplodnaviria","Heunggongvirae","Uroviricota","Caudoviricetes","Lambdavirus","Lambdavirus lambda"]},"proteinDescription":{"recommendedName":{"fullName":{"value":"DNA-packaging protein FI"}}},"genes":[{"geneName":{"value":"Fi"},"orderedLocusNames":[{"value":"lambdap09"}]}],"sequence":{"value":"MTKDELIARLRSLGEQLNRDVSLTGTKEELALRVAELKEELDDTDETAGQDTPLSRENVLTGHENEVGSAQPDTVILDTSELVTVVALVKLHTDALHATRDEPVAFVLPGTAFRVSAGVAAEMTERGLARMQ","length":132,"molWeight":14308,"crc64":"4D7640AC64A97BCA","md5":"D6ACC4387DF6BBF3A5106499FD63FE45"},"extraAttributes":{"uniParcId":"UPI0000138CC0"}}]}
            """
            result = json.loads(result)
            if "results" in result and len(result["results"]) > 0:
                entry = result["results"][0]  # get the first entry

                accession = entry.get("primaryAccession", None)
                gene_name = entry.get("genes", [{}])[0].get("geneName", {}).get("value", None) if entry.get("genes") else None
                locus_tag = entry.get("genes", [{}])[0].get("orderedLocusNames", [{}])[0].get("value", None) if entry.get("genes") else None
                protein_name = entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", None)
                organism_name = entry.get("organism", {}).get("scientificName", None)
                mass = entry.get("sequence", {}).get("molWeight", None)
                sequence = entry.get("sequence", {}).get("value", None)
                return dict(accession=accession,gene_name=gene_name,locus_tag=locus_tag,protein_name=protein_name,organism_name=organism_name,mass=mass,sequence=sequence)

            else:
                print("No results found for the query.")
                return None


    except urllib.error.HTTPError as e:
        print(f"Error: {e.code} - {e.reason}")
    return None


