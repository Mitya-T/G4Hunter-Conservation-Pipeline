from Bio import SeqIO

with open("nipah_everything.fa") as infile, open("nipah_completeGnms.fa", "w") as outfile:
    for record in SeqIO.parse(infile, "fasta"):
        if "complete genome" in record.description:
            SeqIO.write(record, outfile, "fasta")
