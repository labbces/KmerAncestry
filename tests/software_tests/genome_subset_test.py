from Bio import SeqIO

input_file = '/Storage/data1/yuri.dantas/data/SP8032-80_final.fasta'
output_file = '/Storage/data1/yuri.dantas/data/SP803280_subset_2.fasta'

for record in SeqIO.parse(input_file, 'fasta'):
	seq_len = len(record.seq)

	start =  seq_len // 2 + 10**6 //2
	end = start + 10**6
	
	record.seq = record.seq[start:end]
	record.id = record.id + '_1Mbp'
	
	SeqIO.write(record, output_file, 'fasta')
	break

