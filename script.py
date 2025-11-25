from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
import os

input_file = "Benin_gpc.mas"
alignment = AlignIO.read(input_file, "fasta")
alignment_length = alignment.get_alignment_length()
num_sequences = len(alignment)


conserved_indices = []
for i in range(alignment_length):
    column = [record.seq[i] for record in alignment]
    if all(residue == column[0] for residue in column):
        conserved_indices.append(i)


records = []
for record in alignment:
    conserved_seq = ''.join([record.seq[i] for i in conserved_indices])
    records.append(SeqRecord(Seq(conserved_seq),
                             id=record.id,
                             description=""))

alignment_conserved = MultipleSeqAlignment(records)

base_name, ext = os.path.splitext(input_file)
output_file = f"{base_name}_conserved{ext}"

AlignIO.write(alignment_conserved, output_file, "fasta")
print(f"Saved all conserved regions in {output_file}")
