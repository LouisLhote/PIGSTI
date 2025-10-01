# scripts/create_hops_config.py
import sys

original = sys.argv[1]           # e.g. config_hops2.0.txt
new_config = sys.argv[2]         # e.g. config_hops_custom.txt
new_pathogen_list = sys.argv[3]  # e.g. hops_pathogen_list.txt

with open(original) as f:
    lines = f.readlines()

with open(new_config, "w") as f:
    for line in lines:
        if line.strip().startswith("pathToList="):
            f.write(f"pathToList={new_pathogen_list}\n")
        else:
            f.write(line)
