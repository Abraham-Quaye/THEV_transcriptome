import pysam
import csv
from collections import defaultdict

# Set the file paths for the bam and gtf files
bam_files = ["results/hisat2/bulk/subsetTHEV_4hrs.bam", \
"results/hisat2/bulk/subsetTHEV_12hrs.bam", \
"results/hisat2/bulk/subsetTHEV_24hrs.bam", \
"results/hisat2/bulk/subsetTHEV_72hrs.bam"]
gtf_file = "raw_files/annotations/thev_predicted_genes.gtf"

# Iterate over the bam files
for bam_file in bam_files:
    # Initialize a dictionary to store the counts and coordinates for each gene
    gene_counts = defaultdict(lambda: {"count": 0, "start": float("inf"), "end": float("-inf"), "size": 0})
    
    # Read in the bam file
    bam = pysam.AlignmentFile(bam_file, "rb")
    
    # Iterate over the records in the gtf file
    with open(gtf_file, "r") as gtf:
        for line in gtf:
            # Skip comment lines
            if line.startswith("#"):
                continue
            
            # Split the line into fields
            fields = line.strip().split("\t")
            
            # Check if the record is a gene feature
            if fields[2] == "gene":
                # Get the gene ID from the attributes field
                attributes = dict(item.split(" ") for item in fields[8].split("; "))
                gene_name = attributes["Name"].strip("\"")
                
                # Get the chromosome, start, and end positions of the gene
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                
                # Update the start and end coordinates for the gene
                gene_counts[gene_name]["start"] = min(gene_counts[gene_name]["start"], start)
                gene_counts[gene_name]["end"] = max(gene_counts[gene_name]["end"], end)
                
                # Calculate the size of the gene using the blockSizes attribute
                block_sizes = [int(size.strip("\"")) for size in attributes["blockSizes"].split(",")]
                gene_counts[gene_name]["size"] = sum(block_sizes)
                
                # Count the number of reads mapping to the gene
                for read in bam.fetch(chrom, start, end):
                    gene_counts[gene_name]["count"] += 1
    
    # Calculate the total number of reads in the bam file
    total_reads = sum(data["count"] for data in gene_counts.values())
    
    # Calculate the TPM for each gene
    total_tpm = 0
    for data in gene_counts.values():
        count = data["count"]
        size = data["size"]
        tpm = (count / size) * 1e6 / total_reads
        data["tpm"] = tpm
        total_tpm += tpm
        
# Calculate the TPM percentage for each gene
    for data in gene_counts.values():
        tpm = data["tpm"]
        data["tpm_percentage"] = (tpm / total_tpm) * 100
        
    
    # Set the output file path for the csv file
    output_file = f"{bam_file}_pileup.csv"
    
    # Write the counts and coordinates to a csv file
    with open(output_file, "w") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["gene_name", "count", "start", "end", "size", "tpm", "tpm_percentage"])
        for gene_name, data in gene_counts.items():
            read_count = data["count"]
            start = data["start"]
            end = data["end"]
            gene_size = data["size"]
            tpm = data["tpm"]
            tpm_percentage = data["tpm_percentage"]
            writer.writerow([gene_name, read_count, start, end, gene_size, tpm, tpm_percentage])
            
