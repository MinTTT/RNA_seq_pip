# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""

# Built-in/Generic Imports
import os
import sys
# […]

# Libs
import pandas as pd
import numpy as np  # Or any other
# […]

# Own modules
import urllib.parse
import urllib.request

"""
the columns titles for programmatic access: https://www.uniprot.org/help/uniprotkb_column_names
Batch retrieval of entries: https://www.uniprot.org/uploadlists/
database identifier https://www.uniprot.org/help/api_idmapping

"""
url = 'https://www.uniprot.org/uniprot/'
# url = 'https://www.uniprot.org/uploadlists/'

# params = {
# 'from': 'ACC+ID',
# 'to': 'ENSEMBL_ID',
# 'format': 'tab',
# 'query': 'P40925 P40926 O43175 Q9UM73 P97793'
# }
params = {'from': 'ACC', 'to': 'GENENAME', 'format': 'tab',
          'query': 'P0A9H9'}
# params = {'query': 'gene:cheZ AND reviewed:yes AND organism:"Escherichia coli (strain K12) [83333]"',
#           'format': 'tab',
#           'columns': 'id,entry_name,reviewed,protein_names,genes,organism,ec'}


data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
    response = f.read()
print(response.decode('utf-8'))
