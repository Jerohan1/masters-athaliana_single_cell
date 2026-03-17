from BCBio import GFF
from Bio import SeqIO

with open("Col-CC.gff", "w") as out_handle:
    for record in SeqIO.parse("Col-CC.gbff", "genbank"):
        GFF.write([record], out_handle)
