# This code was written with the help of ChatGPT

import json
import os
from pathlib import Path

# Define the base directory containing all genus folders
base_dir = 'bateria'

# Output file to store all concatenated FASTA sequences
output_file = 'combined_bacteria_genomes.fna'

# Initialize a counter for the number of sequences combined
sequence_count = 0

# Open the output file for writing
with open(output_file, 'w') as outfile:
    # Recursively search for `dataset_catalog.json` in each genus directory
    for genus_dir in Path(base_dir).rglob('data/dataset_catalog.json'):
        # Load each `dataset_catalog.json`
        with open(genus_dir, 'r') as f:
            catalog = json.load(f)
        
        # Find the directory containing the genome files
        genus_data_dir = genus_dir.parent
        
        # Iterate over each assembly listed in the catalog
        for assembly in catalog['assemblies']:
            for file_info in assembly['files']:
                # Check if the file is a genomic FASTA file
                if file_info['fileType'] == 'GENOMIC_NUCLEOTIDE_FASTA':
                    # Construct the full path to the .fna file
                    fasta_path = genus_data_dir / file_info['filePath']
                    
                    # Open and read each FASTA file, then append it to the output file
                    with open(fasta_path, 'r') as fasta_file:
                        # Read the file content
                        fasta_content = fasta_file.read()
                        
                        # Write the content to the output file
                        outfile.write(fasta_content)
                        
                        # Increment the sequence count based on the number of '>' symbols (FASTA headers)
                        sequence_count += fasta_content.count('>')

# Output the number of sequences combined
print(f"All FASTA files have been concatenated into {output_file}.")
print(f"Total number of sequences combined: {sequence_count}")