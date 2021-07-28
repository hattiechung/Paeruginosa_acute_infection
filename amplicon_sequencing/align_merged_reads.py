import os, subprocess
import pandas as pd

def get_path():
	if os.path.exists('/Volumes'):
		pathprefix = '/Volumes/ahg_regevdata'
	else:
		pathprefix = '/ahg/regevdata'
		
	filepath = '%s/projects/olfaction_sc/PA_amplicon/fastq/PA_amplicon' %pathprefix
	tablespath = '%s/projects/olfaction_sc/PA_amplicon' %pathprefix
	return filepath, tablespath

def load_sample_tables(): 
	_, tablespath = get_path()
	df_min_len = pd.read_csv('%s/gene_min_length.csv' %tablespath)
	df_primers = pd.read_csv('%s/gene_primers.csv' %tablespath) 
	return df_min_len, df_primers

def load_colony_freqs(): 
	_, tablespath = get_path()
	df_colony_freqs = pd.read_csv('%s/isolate_mut_freqs.csv' %tablespath)
	return df_colony_freqs

def get_sample_sheet():
	filepath, tablespath = get_path() 
	df_sample = pd.read_csv('%s/sample_sheet.csv' %tablespath)
	fastqs = os.listdir(filepath)
	fastqs.sort()
	for f in fastqs:
		if not f.startswith('blank'):
			locus, sample, spacer, readside, extension = f.split('_')
			sample_row = (df_sample['literal_gene']==locus) & (df_sample['sample']==sample)
			if readside=='R1':
				df_sample.loc[sample_row,'read1'] = f
			elif readside=='R2': 
				df_sample.loc[sample_row,'read2'] = f
	return df_sample

def get_length_thresh(row): 
    a = df_min_len[df_min_len['Gene']==row['gene']]['Min_len_R1']
    b = df_min_len[df_min_len['Gene']==row['gene']]['Min_len_R2']
    return int(a), int(b)

def bowtie2_align_merged_read(sdf): 
	filepath, tablespath = get_path() 
	bowtie2_ref = '/ahg/regevdata/projects/olfaction_sc/PA_amplicon/amplicon_ref'
	
	for idx, row in sdf.iterrows(): 
		dir_name = row['gene']+'_'+row['sample_type']
		out_dir = tablespath+'/cutadapt/'+dir_name
		print(out_dir)
		if not os.path.isfile(out_dir+'/merged.sam'): 
			p1 = subprocess.run(['bowtie2', '--local', '-x', bowtie2_ref, 
							'-U', out_dir+'/merged.fastq', 
							'-S', out_dir+'/merged.sam'])

def filter_mapped(sdf): 
	filepath, tablespath = get_path() 
	
	bowtie2_ref = '/ahg/regevdata/projects/olfaction_sc/PA_amplicon/amplicon_ref'
	fasta_ref = '/ahg/regevdata/projects/olfaction_sc/PA_amplicon/targets.fasta'
	for idx, row in sdf.iterrows(): 
		dir_name = row['gene']+'_'+row['sample_type']
		out_dir = tablespath+'/cutadapt/'+dir_name
		print(out_dir)
		if os.path.isfile(out_dir+'/merged.sam'): 
# 			p1 = subprocess.run(['samtools', 'view', 
# 							out_dir+'/merged.sam'
# 							'-o', out_dir+'/merged_filtered.bam'], stdout=subprocess.PIPE)
	
			print('sorting SAM')
			p2 = subprocess.run(['samtools', 'sort', 
							'-o', out_dir+'/sorted.bam',
							out_dir+'/merged.sam'],
							stdout=subprocess.PIPE)

			print('counting UMIs')
			p3 = subprocess.run(['cut', '-f', '10', out_dir+'/merged.sam', '\|',
								'cut', '-c1-8', '\|', 
								'sort', '-r', '\|', 'uniq', '-c', '\|', 'sort', '-nrk1,1', '>>', 
								out_dir+'/UMIs.txt'],
								stdout=subprocess.PIPE)
			
			print('mpileup')
			p4 = subprocess.run(['samtools', 'mpileup', '-O', '-s'
							'-f', fasta_ref, 
							'-o', out_dir+'/sorted.pileup',
							out_dir+'/sorted.bam'],
							stdout=subprocess.PIPE)

			p5 = subprocess.run(['bcftools', 'view', '-bvcg',
							out_dir+'/sorted.bcf', '>>', 
							out_dir+'/variants.bcf'],
							shell=True)
			# p6 = subprocess.run(['bcftools', 'view', 
# 							out_dir+'/variants.bcf', '>>', 
# 							out_dir+'/variants_filtered.vcf'],
# 							shell=True)

def count_UMIs(sdf): 
	
	filepath, tablespath = get_path() 
	df_min_len, primer_table = load_sample_tables()
	
	for idx, row in sdf.iterrows(): 
		dir_name = row['gene']+'_'+row['sample_type']
		out_dir = tablespath+'/cutadapt/'+dir_name

		if os.path.isfile(out_dir+'/merged.sam'): 
			l_R1, _ = get_length_thresh(row)
			grab_len = l_R1+6
			
			# ___ Extract aligned reads from SAM file ___ #
			
			outfile = out_dir+'/mapped_reads.txt'
			if not os.path.isfile(outfile):
				print('Extracting mapped reads')
				f = open(outfile, 'w')
				cmds = ['grep', '-n', '"'+row['gene_grep']+'"', out_dir+'/merged.sam', '|', 'cut', '-f', '10']
				subprocess.Popen(' '.join(cmds), shell=True, stdout=f)
				f.close()
			
			# ___ Extract UMIs ___ #
			
			f_UMI1 = out_dir+'/UMI1.txt'
			if not os.path.isfile(f_UMI1): 
				print('Extracting UMI 1')
				f = open(f_UMI1, 'w')
				cmds = ['cat', out_dir+'/mapped_reads.txt', '|', 'cut', '-c1-8']
				subprocess.Popen(' '.join(cmds), shell=True, stdout=f)
				f.close()
					
			f_UMI2 = out_dir+'/UMI2.txt'	
			if not os.path.isfile(f_UMI2):
				print('Extracting UMI 2')
				f2 = open(f_UMI2, 'w')
				cmds = ['cat', out_dir+'/mapped_reads.txt', '|', 'rev', '|', 'cut', '-c1-8']
				subprocess.Popen(' '.join(cmds), shell=True, stdout=f2)
				f2.close()
					
			f_UMI_c = out_dir+'/UMI_combined.txt'
			f_UMI_s = out_dir+'/UMIstats.txt'
			if not os.path.isfile(f_UMI_s):
				print('Combining joint UMIs')
				f3 = open(f_UMI_c, 'w')
				
				# Merge UMIs
				UMI1 = open(f_UMI1,'r').readlines()
				UMI2 = open(f_UMI2,'r').readlines()
				for umi1, umi2 in zip(UMI1, UMI2): 
					combined_UMI = umi1.strip()+umi2.strip()
					f3.write('%s\n' %combined_UMI)
				f3.close()
				
				f4 = open(f_UMI_s, 'w')
				cmds = ['cat', out_dir+'/UMI_combined.txt', '|', 'sort', '-r', '|', 
						'uniq', '-c', '|', 'sort', '-nrk1,1']
				subprocess.Popen(' '.join(cmds), shell=True, stdout=f4)
				f4.close()

def get_unique_nonUMI_reads(sdf): 
	filepath, tablespath = get_path() 
	df_min_len, primer_table = load_sample_tables()
	
	for idx, row in sdf.iterrows(): 
		dir_name = row['gene']+'_'+row['sample_type']
		out_dir = tablespath+'/cutadapt/'+dir_name

		if os.path.isfile(out_dir+'/merged.sam'): 

			cmds = ['grep', '-n', row['gene'], out_dir+'/merged.sam', '|',
								'cut', '-f', '10', '|',
								'cut', '-c9-'+str(grab_len), '|', 
								'sort', '-r', '|', 'uniq', '-c', '|', 'sort', '-nrk1,1', '>', 
								out_dir+'/unique_reads.txt']
			print(' '.join(cmds))
			p3 = subprocess.run(cmds, stdout=subprocess.PIPE)

def extract_variant_freqs(sdf): 
	# return tables (df) per gene/locus, with merged frequencies across isogenic controls and t1,t2
	import regex
	filepath, tablespath = get_path() 
	df_min_len, primer_table = load_sample_tables()

	genes = set(sdf['gene'])
	GENE_CTS = {gene: pd.DataFrame() for gene in genes}
	
	for idx, row in sdf.iterrows(): 
	
		dir_name = row['gene']+'_'+row['sample_type']
		out_dir = tablespath+'/cutadapt/'+dir_name
		
		var_start = primer_table[primer_table['Gene']==row['gene']]['var_start'].values[0]
		var_end = primer_table[primer_table['Gene']==row['gene']]['var_end'].values[0]
		
		gene_FW = primer_table[primer_table['Gene']==row['gene']]['FW'].values[0]
		
		UMIs = open(out_dir+'/UMI_combined.txt', 'r').readlines()
		reads = open(out_dir+'/mapped_reads.txt', 'r').readlines()
		
		# make variant-UMI dictionary
		var_dict = {}
		for umi, seq in zip(UMIs, reads):
			m = regex.search("(%s){s<=2}" %(gene_FW), seq)		# s<=2 vs e<=2 substitution vs. elastic
			if m: 
				m_start = m.span()[0]
				# get variants using coords relative to FW primer start 
				var_seq = seq[m_start+var_start:m_start+var_end+1]
					
				if var_seq not in var_dict.keys(): 
					var_dict[var_seq] = [umi.strip()]
				else: 
					var_dict[var_seq].append(umi.strip())
		
		# count UMI freqs per variant
		var_ct = {}
		for k, vals in var_dict.items(): 
			var_ct[k] = len(set(vals))
		cts = pd.Series(var_ct, name=row['sample_type'], index=var_ct.keys())
		
		if out_dir=='ampR_MUT': print(cts)
		
		# add to master table	
# 		GENE_CTS[row['gene']][row['sample_type']] = cts
		GENE_CTS[row['gene']] = pd.concat([GENE_CTS[row['gene']], cts.to_frame()], axis=1, join='outer')
	
	# sort and filter
	for gene in genes: 
		GENE_CTS[gene].sort_values(["t1", "t2"], ascending=(False, False), inplace=True)
		GENE_CTS[gene].drop(GENE_CTS[gene][GENE_CTS[gene].sum(axis=1)<5].index, inplace = True) 
		
	return GENE_CTS
		
if __name__ == '__main__': 
	samples_df = get_sample_sheet()
	
# 	bowtie2_align_merged_read(samples_df)
# 	count_UMIs(samples_df)
	extract_variant_freqs(samples_df)
