# Copyright © 2022, Shenghao Cao, Mengtian Li & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China 

import sys
import re

with open(sys.argv[1]) as fin:
	with open(sys.argv[2],'w+') as fout:
		
		for line in fin:
			
			record = line.strip()
			string = record.split('|')
			num = len(string)
			
			if line.startswith('>'):
				
				scaff_num = re.match(r'>scaffold(\d+)',string[0]).group(1)
				scaff_size = re.match(r'size(\d+)',string[1]).group(1)
				scaff_tigs = re.match(r'tigs(\d+)',string[2]).group(1)
				print("~~~~~",file=fout)
				print(scaff_num+"\t"+scaff_size+"\t"+scaff_tigs,file=fout)
			
			elif num > 1:

				links = '0'
				gaps = '0'
				gaps_lowbound = '0'
				gaps_upbound = '0'
				merged = '0'

				orien = re.match(r'(\w)_tig(\d+)',string[0]).group(1)
				contig_num = re.match(r'(\w)_tig(\d+)',string[0]).group(2)
				contig_size = re.match(r'size(\d+)',string[1]).group(1)

				if num >= 4:

					links = re.match(r'links(\d+)',string[2]).group(1)
					gaps = re.match(r'gaps(-*\d+)',string[3]).group(1)
					gaps_lowbound = gaps
					gaps_upbound = gaps
					if num >= 5:
						merged = re.match(r'merged(\d+)',string[4]).group(1)

				print(contig_num+'\t'+orien+'\t'+contig_size+'\t'+links+'\t'+gaps+'\t'+gaps_lowbound+'\t'+gaps_upbound+'\t'+merged,file=fout)