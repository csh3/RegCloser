# Run RegCloser
python ~/RegCloser/RunPipeline.py -p prerequisite -g draft_genome.fasta -d output -t 48

# Run Phrap for gap closing
python ~/RegCloser/RunPipeline_phrap.py -p prerequisite -g draft_genome.fasta -d output_phrap -t 48
