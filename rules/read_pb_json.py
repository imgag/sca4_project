#! python 
# Read and collate json samples
import json
import glob
import os
from pathlib import Path

json_files = glob.glob('**/datastore*.json', recursive=True)

j = {}
for json_path in json_files:
    json_full_path = (Path(json_path).resolve())
    with open(json_full_path, 'r') as json_file:
        j[str(json_full_path)] = json.loads(json_file.read())

print(json.dumps(j, indent = 4, sort_keys = True))