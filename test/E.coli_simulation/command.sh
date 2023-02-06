# Run RegCloser
python ~/RegCloser/RunPipeline.py -p prerequisite -g draft_genome.fasta -d iter-1 -w -t 48
python ~/RegCloser/RunPipeline.py -p prerequisite -g iter-1/output_genome.fasta -d iter-2 -w -t 48

# Run Phrap for gap closing
python ~/RegCloser/RunPipeline_phrap.py -p prerequisite -g draft_genome.fasta -d iter-1_phrap -t 48
python ~/RegCloser/RunPipeline_phrap.py -p prerequisite -g iter-1_phrap/output_genome.fasta -d iter-2_phrap -t 48
