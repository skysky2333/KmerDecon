# This code was written with the help of ChatGPT

import subprocess

# Define all taxons
all_taxons = [
    "Arthrobacter",
    "Burkholderia",
    "Chryseobacterium",
    "Ochrobactrum",
    "Pseudomonas",
    "Ralstonia",
    "1827", #Rhodococcus
    "Sphingomonas",
    "Corynebacterium",
    "Propionibacterium",
    "Streptococcus"
]

# Iterate through taxons and execute the command
for taxon in all_taxons:

    # Construct the command
    command = [
        "./datasets",
        "download",
        "genome",
        "taxon",
        taxon,
        "--reference",
        "--assembly-level",
        "complete,chromosome",
        "--filename",
        taxon + ".zip"
    ]

    print(f"Running command for taxon: {taxon}")
    
    # Execute the command
    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error occurred while processing {taxon}: {e}")