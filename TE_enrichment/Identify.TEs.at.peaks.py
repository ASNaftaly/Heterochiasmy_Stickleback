#pull observed counts of peaks at TEs broad class analysis first
#uses logic from random permutations files
#to run script: python3 Identify.TEs.at.peaks.py <peaks xls file> <TEs .out file> <output: DNA TEs at peaks> <output: LINE TEs at peaks> <output: SINE TEs at peaks> <output: LTR TEs at peaks> <output: RC TEs at peaks>
#Author: Alice Naftaly, Sept 2020

import sys

#pulls ATAC seq peaks from file (xls) from Macs2
#returns dictionary with following format:
#key = chr number
#dict value = [peak start, peak end]
def pull_peaks():
	peaks_file = sys.argv[1]
	peaks_dict = {}
	with open(peaks_file, 'r') as peaks:
		for line in peaks:
			if line.startswith("chr"):
				new_line = line.split()
				if new_line[0] == "chr":
					continue
				else:
					chr_num = new_line[0]
					start = new_line[1]
					end = new_line[2]
					dict_value = [start,end]
					if chr_num in peaks_dict:
						peaks_dict[chr_num].append(dict_value)
					elif chr_num not in peaks_dict:
						peaks_dict.update({chr_num:[dict_value]})
	return peaks_dict


#pulls TEs from repeat masker file
#format of repeat masker file: SW.repeat.score  Perc.Div    Perc.Del    Perc.Ins    Query.Seq   Query.Start.pos Query.End.pos   Query (left)    +/C     Matching.repeat Repeat.Class/Family     Repeat.Start.pos    Repeat.End.pos  Repeat(left)    ID
#returns dictionary with following format:
#key = feature (DNA, LINES, SINES, LTRs, RC)
#dict values = [chr number, te name, te start pos, te end pos]
def pull_TEs():
	te_file = sys.argv[2]
	te_dict = {}
	with open(te_file, 'r') as te:
		for line in te:
			if len(line) != 124 and len(line) != 131 and len(line) != 1:
				new_line = line.split()
				chr_num = new_line[4]
				feature = new_line[10]
				if feature.startswith("DNA"):
					final_feature = "DNA"
					start = new_line[5]
					end = new_line[6]
					te_name = new_line[10]
					dict_value = [chr_num, te_name,start, end]
					if final_feature in te_dict:
						te_dict[final_feature].append(dict_value)
					elif final_feature not in te_dict:
						te_dict.update({final_feature:[dict_value]})
				elif feature.startswith("LINE"):
					final_feature = "LINES"
					start = new_line[5]
					end = new_line[6]
					te_name = new_line[10]
					dict_value = [chr_num, te_name,start, end]
					if final_feature in te_dict:
						te_dict[final_feature].append(dict_value)
					elif final_feature not in te_dict:
						te_dict.update({final_feature:[dict_value]})
				elif feature.startswith("SINE"):
					final_feature = "SINES"
					start = new_line[5]
					end = new_line[6]
					te_name = new_line[10]
					dict_value = [chr_num, te_name,start, end]
					if final_feature in te_dict:
						te_dict[final_feature].append(dict_value)
					elif final_feature not in te_dict:
						te_dict.update({final_feature:[dict_value]})
				elif feature.startswith("LTR"):
					final_feature = "LTRs"
					start = new_line[5]
					end = new_line[6]
					te_name = new_line[10]
					dict_value = [chr_num, te_name,start, end]
					if final_feature in te_dict:
						te_dict[final_feature].append(dict_value)
					elif final_feature not in te_dict:
						te_dict.update({final_feature:[dict_value]})
				elif feature.startswith("RC"):
					final_feature = "RC"
					start = new_line[5]
					end = new_line[6]
					te_name = new_line[10]
					dict_value = [chr_num, te_name,start, end]
					if final_feature in te_dict:
						te_dict[final_feature].append(dict_value)
					elif final_feature not in te_dict:
						te_dict.update({final_feature:[dict_value]})
	return te_dict


#Compare peaks with each category of TEs
#for each function, the dictionary returned is:
#key = type of te/peak overlap (tops, tip, pit, tope)
#dictionary value = chr_num, te_name, te_start, te_end, peak_start, peak_end
def compare_DNA():
	peaks = pull_peaks()
	TEs = pull_TEs()
	final_te_dict = {}
	for value in TEs["DNA"]:
		chr_num = value[0]
		te_name = value[1]
		te_start_final = value[2]
		te_end_final = value[3]
		final = [te_name, te_start_final,te_end_final]
		if chr_num in final_te_dict:
			final_te_dict[chr_num].append(final)
		elif chr_num not in final_te_dict:
			final_te_dict.update({chr_num:[final]})
	DNA_TE_dict = {}
	for key in peaks:
		if key in final_te_dict:
			peaks_value = peaks[key]
			te_values = final_te_dict[key]
			chr_num = key
			for val in te_values:
				te_name = str(val[0])
				te_start = int(val[1])
				te_end = int(val[2])
				for v in peaks_value:
					peak_start = int(v[0])
					peak_end = int(v[1])
					#pit = peak in TE
					if te_start <= peak_start and te_start <= peak_end and te_end >= peak_start and te_end >= peak_end:
						dict_value = ["pit", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in DNA_TE_dict:
							DNA_TE_dict[chr_num].append(dict_value)
						elif chr_num not in DNA_TE_dict:
							DNA_TE_dict.update({chr_num:[dict_value]})
					#pots = peak overlaps TE start
					elif te_start > peak_start and te_start < peak_end and te_end >= peak_start and te_end >= peak_end:
						dict_value = ["pots", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in DNA_TE_dict:
							DNA_TE_dict[chr_num].append(dict_value)
						elif chr_num not in DNA_TE_dict:
							DNA_TE_dict.update({chr_num:[dict_value]})
					#pote = peak overlaps TE end
					elif te_start <= peak_start and te_start <= peak_end and te_end > peak_start and te_end < peak_end:
						dict_value = ["pote", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in DNA_TE_dict:
							DNA_TE_dict[chr_num].append(dict_value)
						elif chr_num not in DNA_TE_dict:
							DNA_TE_dict.update({chr_num:[dict_value]})
					#pot = peak overlaps TE
					elif te_start <= peak_start and te_start <= peak_end and te_end > peak_start and te_end < peak_end:
						dict_value = ["pot", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in DNA_TE_dict:
							DNA_TE_dict[chr_num].append(dict_value)
						elif chr_num not in DNA_TE_dict:
							DNA_TE_dict.update({chr_num:[dict_value]})
	return DNA_TE_dict

#LINES
#for each function, the dictionary returned is:
#key = type of te/peak overlap (tops, tip, pit, tope)
#dictionary value = chr_num, te_name, te_start, te_end, peak_start, peak_end
def compare_LINES():
	peaks = pull_peaks()
	TEs = pull_TEs()
	final_te_dict = {}
	for value in TEs["LINES"]:
		chr_num = value[0]
		te_name = value[1]
		te_start_final = value[2]
		te_end_final = value[3]
		final = [te_name, te_start_final,te_end_final]
		if chr_num in final_te_dict:
			final_te_dict[chr_num].append(final)
		elif chr_num not in final_te_dict:
			final_te_dict.update({chr_num:[final]})
	LINE_TE_dict = {}
	for key in peaks:
		if key in final_te_dict:
			peaks_value = peaks[key]
			te_values = final_te_dict[key]
			chr_num = key
			for val in te_values:
				te_name = str(val[0])
				te_start = int(val[1])
				te_end = int(val[2])
				for v in peaks_value:
					peak_start = int(v[0])
					peak_end = int(v[1])
					#pit = peak in TE
					if te_start <= peak_start and te_start <= peak_end and te_end >= peak_start and te_end >= peak_end:
						dict_value = ["pit", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in LINE_TE_dict:
							LINE_TE_dict[chr_num].append(dict_value)
						elif chr_num not in LINE_TE_dict:
							LINE_TE_dict.update({chr_num:[dict_value]})
					#pots = peak overlaps TE start
					elif te_start > peak_start and te_start < peak_end and te_end >= peak_start and te_end >= peak_end:
						dict_value = ["pots", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in LINE_TE_dict:
							LINE_TE_dict[chr_num].append(dict_value)
						elif chr_num not in LINE_TE_dict:
							LINE_TE_dict.update({chr_num:[dict_value]})
					#pote = peak overlaps TE end
					elif te_start <= peak_start and te_start <= peak_end and te_end > peak_start and te_end < peak_end:
						dict_value = ["pote", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in LINE_TE_dict:
							LINE_TE_dict[chr_num].append(dict_value)
						elif chr_num not in LINE_TE_dict:
							LINE_TE_dict.update({chr_num:[dict_value]})
					#pot = peak overlaps TE
					elif te_start <= peak_start and te_start <= peak_end and te_end > peak_start and te_end < peak_end:
						dict_value = ["pot", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in LINE_TE_dict:
							LINE_TE_dict[chr_num].append(dict_value)
						elif chr_num not in LINE_TE_dict:
							LINE_TE_dict.update({chr_num:[dict_value]})
	return LINE_TE_dict

#SINES
#for each function, the dictionary returned is:
#key = type of te/peak overlap (tops, tip, pit, tope)
#dictionary value = chr_num, te_name, te_start, te_end, peak_start, peak_end
def compare_SINES():
	peaks = pull_peaks()
	TEs = pull_TEs()
	final_te_dict = {}
	for value in TEs["SINES"]:
		chr_num = value[0]
		te_name = value[1]
		te_start_final = value[2]
		te_end_final = value[3]
		final = [te_name, te_start_final,te_end_final]
		if chr_num in final_te_dict:
			final_te_dict[chr_num].append(final)
		elif chr_num not in final_te_dict:
			final_te_dict.update({chr_num:[final]})
	SINE_TE_dict = {}
	for key in peaks:
		if key in final_te_dict:
			peaks_value = peaks[key]
			te_values = final_te_dict[key]
			chr_num = key
			for val in te_values:
				te_name = str(val[0])
				te_start = int(val[1])
				te_end = int(val[2])
				for v in peaks_value:
					peak_start = int(v[0])
					peak_end = int(v[1])
					#pit = peak in TE
					if te_start <= peak_start and te_start <= peak_end and te_end >= peak_start and te_end >= peak_end:
						dict_value = ["pit", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in SINE_TE_dict:
							SINE_TE_dict[chr_num].append(dict_value)
						elif chr_num not in SINE_TE_dict:
							SINE_TE_dict.update({chr_num:[dict_value]})
					#pots = peak overlaps TE start
					elif te_start > peak_start and te_start < peak_end and te_end >= peak_start and te_end >= peak_end:
						dict_value = ["pots", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in SINE_TE_dict:
							SINE_TE_dict[chr_num].append(dict_value)
						elif chr_num not in SINE_TE_dict:
							SINE_TE_dict.update({chr_num:[dict_value]})
					#pote = peak overlaps TE end
					elif te_start <= peak_start and te_start <= peak_end and te_end > peak_start and te_end < peak_end:
						dict_value = ["pote", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in SINE_TE_dict:
							SINE_TE_dict[chr_num].append(dict_value)
						elif chr_num not in SINE_TE_dict:
							SINE_TE_dict.update({chr_num:[dict_value]})
					#pot = peak overlaps TE
					elif te_start <= peak_start and te_start <= peak_end and te_end > peak_start and te_end < peak_end:
						dict_value = ["pot", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in SINE_TE_dict:
							SINE_TE_dict[chr_num].append(dict_value)
						elif chr_num not in SINE_TE_dict:
							SINE_TE_dict.update({chr_num:[dict_value]})
	return SINE_TE_dict


#LTRs
#for each function, the dictionary returned is:
#key = type of te/peak overlap (tops, tip, pit, tope)
#dictionary value = chr_num, te_name, te_start, te_end, peak_start, peak_end
def compare_LTRs():
	peaks = pull_peaks()
	TEs = pull_TEs()
	final_te_dict = {}
	for value in TEs["LTRs"]:
		chr_num = value[0]
		te_name = value[1]
		te_start_final = value[2]
		te_end_final = value[3]
		final = [te_name, te_start_final,te_end_final]
		if chr_num in final_te_dict:
			final_te_dict[chr_num].append(final)
		elif chr_num not in final_te_dict:
			final_te_dict.update({chr_num:[final]})
	LTR_TE_dict = {}
	for key in peaks:
		if key in final_te_dict:
			peaks_value = peaks[key]
			te_values = final_te_dict[key]
			chr_num = key
			for val in te_values:
				te_name = str(val[0])
				te_start = int(val[1])
				te_end = int(val[2])
				for v in peaks_value:
					peak_start = int(v[0])
					peak_end = int(v[1])
					#pit = peak in TE
					if te_start <= peak_start and te_start <= peak_end and te_end >= peak_start and te_end >= peak_end:
						dict_value = ["pit", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in LTR_TE_dict:
							LTR_TE_dict[chr_num].append(dict_value)
						elif chr_num not in LTR_TE_dict:
							LTR_TE_dict.update({chr_num:[dict_value]})
					#pots = peak overlaps TE start
					elif te_start > peak_start and te_start < peak_end and te_end >= peak_start and te_end >= peak_end:
						dict_value = ["pots", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in LTR_TE_dict:
							LTR_TE_dict[chr_num].append(dict_value)
						elif chr_num not in LTR_TE_dict:
							LTR_TE_dict.update({chr_num:[dict_value]})
					#pote = peak overlaps TE end
					elif te_start <= peak_start and te_start <= peak_end and te_end > peak_start and te_end < peak_end:
						dict_value = ["pote", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in LTR_TE_dict:
							LTR_TE_dict[chr_num].append(dict_value)
						elif chr_num not in LTR_TE_dict:
							LTR_TE_dict.update({chr_num:[dict_value]})
					#pot = peak overlaps TE
					elif te_start <= peak_start and te_start <= peak_end and te_end > peak_start and te_end < peak_end:
						dict_value = ["pot", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in LTR_TE_dict:
							LTR_TE_dict[chr_num].append(dict_value)
						elif chr_num not in LTR_TE_dict:
							LTR_TE_dict.update({chr_num:[dict_value]})
	return LTR_TE_dict

#RC
#for each function, the dictionary returned is:
#key = type of te/peak overlap (tops, tip, pit, tope)
#dictionary value = chr_num, te_name, te_start, te_end, peak_start, peak_end
def compare_RC():
	peaks = pull_peaks()
	TEs = pull_TEs()
	final_te_dict = {}
	for value in TEs["RC"]:
		chr_num = value[0]
		te_name = value[1]
		te_start_final = value[2]
		te_end_final = value[3]
		final = [te_name, te_start_final,te_end_final]
		if chr_num in final_te_dict:
			final_te_dict[chr_num].append(final)
		elif chr_num not in final_te_dict:
			final_te_dict.update({chr_num:[final]})
	RC_TE_dict = {}
	for key in peaks:
		if key in final_te_dict:
			peaks_value = peaks[key]
			te_values = final_te_dict[key]
			chr_num = key
			for val in te_values:
				te_name = str(val[0])
				te_start = int(val[1])
				te_end = int(val[2])
				for v in peaks_value:
					peak_start = int(v[0])
					peak_end = int(v[1])
					#pit = peak in TE
					if te_start <= peak_start and te_start <= peak_end and te_end >= peak_start and te_end >= peak_end:
						dict_value = ["pit", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in RC_TE_dict:
							RC_TE_dict[chr_num].append(dict_value)
						elif chr_num not in RC_TE_dict:
							RC_TE_dict.update({chr_num:[dict_value]})
					#pots = peak overlaps TE start
					elif te_start > peak_start and te_start < peak_end and te_end >= peak_start and te_end >= peak_end:
						dict_value = ["pots", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in RC_TE_dict:
							RC_TE_dict[chr_num].append(dict_value)
						elif chr_num not in RC_TE_dict:
							RC_TE_dict.update({chr_num:[dict_value]})
					#pote = peak overlaps TE end
					elif te_start <= peak_start and te_start <= peak_end and te_end > peak_start and te_end < peak_end:
						dict_value = ["pote", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in RC_TE_dict:
							RC_TE_dict[chr_num].append(dict_value)
						elif chr_num not in RC_TE_dict:
							RC_TE_dict.update({chr_num:[dict_value]})
					#pot = peak overlaps TE
					elif te_start <= peak_start and te_start <= peak_end and te_end > peak_start and te_end < peak_end:
						dict_value = ["pot", te_name, te_start, te_end, peak_start, peak_end]
						if chr_num in RC_TE_dict:
							RC_TE_dict[chr_num].append(dict_value)
						elif chr_num not in RC_TE_dict:
							RC_TE_dict.update({chr_num:[dict_value]})
	return RC_TE_dict

#write output files that show chr number, TE name, te start, te end, peak start peak end for each broad class of TE families
#will have separate output for each broad TE class
def write_dna():
	dna_tes = compare_DNA()
	chrs = ["chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXIX","chrXX","chrXXI"]
	DNA_output = sys.argv[3]
	with open(DNA_output, 'a') as DNA_out:
		#total of 7 columns for main part of document
		header_DNA = "DNA.TEs overlapping sigpeaks\nChr\tOverlap.Type\tTE.name\tTE.start\tTE.end\tPeak.Start\tPeak.end\n"
		DNA_out.write(header_DNA)
		for value in chrs:
			if value in dna_tes:
				dna_single_chr = dna_tes[value]
				for val in dna_single_chr:
					dna_final = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(value),str(val[0]),str(val[1]),str(val[2]),str(val[3]),str(val[4]),str(val[5]))
					DNA_out.write(dna_final)
			else:
				continue

def write_lines():
	lines_tes = compare_LINES()
	chrs = ["chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXIX","chrXX","chrXXI"]
	lines_output = sys.argv[4]
	with open(lines_output, 'a') as LINES_out:
		#total of 7 columns for main part of document
		header_LINES = "LINES.TEs overlapping sigpeaks\nChr\tOverlap.Type\tTE.name\tTE.start\tTE.end\tPeak.Start\tPeak.end\n"
		LINES_out.write(header_LINES)
		for value in chrs:
			if value in lines_tes:
				lines_single_chr = lines_tes[value]
				for val in lines_single_chr:
					lines_final = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(value),str(val[0]),str(val[1]),str(val[2]),str(val[3]),str(val[4]),str(val[5]))
					LINES_out.write(lines_final)
			else:
				continue

def write_sines():
	sines_tes = compare_SINES()
	chrs = ["chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXIX","chrXX","chrXXI"]
	sines_output = sys.argv[5]
	with open(sines_output, 'a') as SINES_out:
		#total of 7 columns for main part of document
		header_SINES = "SINES.TEs overlapping sigpeaks\nChr\tOverlap.Type\tTE.name\tTE.start\tTE.end\tPeak.Start\tPeak.end\n"
		SINES_out.write(header_SINES)
		for value in chrs:
			if value in sines_tes:
				sines_single_chr = sines_tes[value]
				for val in sines_single_chr:
					sines_final = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(value),str(val[0]),str(val[1]),str(val[2]),str(val[3]),str(val[4]),str(val[5]))
					SINES_out.write(sines_final)
			else:
				continue

def write_ltrs():
	ltrs_tes = compare_LTRs()
	chrs = ["chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXIX","chrXX","chrXXI"]
	ltrs_output = sys.argv[6]
	with open(ltrs_output, 'a') as LTRs_out:
		#total of 7 columns for main part of document
		header_LTRs = "LTRs.TEs overlapping sigpeaks\nChr\tOverlap.Type\tTE.name\tTE.start\tTE.end\tPeak.Start\tPeak.end\n"
		LTRs_out.write(header_LTRs)
		for value in chrs:
			if value in ltrs_tes:
				ltrs_single_chr = ltrs_tes[value]
				for val in ltrs_single_chr:
					ltrs_final = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(value),str(val[0]),str(val[1]),str(val[2]),str(val[3]),str(val[4]),str(val[5]))
					LTRs_out.write(ltrs_final)
			else:
				continue

def write_rc():
	rc_tes = compare_RC()
	chrs = ["chrI","chrII","chrIII","chrIV","chrV","chrVI","chrVII","chrVIII","chrIX","chrX","chrXI","chrXII","chrXIII","chrXIV","chrXV","chrXVI","chrXVII","chrXVIII","chrXIX","chrXX","chrXXI"]
	rc_output = sys.argv[7]
	with open(rc_output, 'a') as RC_out:
		#total of 7 columns for main part of document
		header_RC = "RC.TEs overlapping sigpeaks\nChr\tOverlap.Type\tTE.name\tTE.start\tTE.end\tPeak.Start\tPeak.end\n"
		RC_out.write(header_RC)
		for value in chrs:
			if value in rc_tes:
				rc_single_chr = rc_tes[value]
				for val in rc_single_chr:
					rc_final = "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(value),str(val[0]),str(val[1]),str(val[2]),str(val[3]),str(val[4]),str(val[5]))
					RC_out.write(rc_final)
			else:
				continue


#write alloutput files
def write_all():
	dna = write_dna()
	lines = write_lines()
	sines = write_sines()
	ltrs = write_ltrs()
	rc = write_rc()

write_all()
