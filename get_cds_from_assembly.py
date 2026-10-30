### V0.3 ###
### Shakunthala Natarajan ###
"""
feature requests and bug reports to s64snata@uni-bonn.de
"""

__version__='0.3'

__usage__ = """
			python3 get_cds_from_assembly.py
			--assembly <full path to file or folder of assembly files>
			--gff <full path to file or folder of gff files>
			--out <full path to output directory>

			NOTE_1: Assembly files should have extension of .fa, .fna or .fasta (optionally followed by .gz)
			NOTE_2: GFF files should have extension of .gff or .gff3 (optionally followed by .gz)
			NOTE_3: An assembly file and a GFF file are paired when their names are identical after removing the extension (case ignored),
					e.g. species1.fasta.gz and species1.gff3
			NOTE_4: GFF files without CDS features are skipped with a warning, because exon-based sequences would include UTRs
			NOTE_5: GFF seqids (column 1) must match the FASTA headers (first word only). Transcripts on contigs that are missing
					in the FASTA are skipped; if no seqid matches at all, the species is skipped with an error message

			optional:
					--gff_config <full path to txt file containing the gff parameters to be used separated by spaces or tabs;
					Needs 4 columns in the order 
					(i) base file name - same as the base name you use for the gff file (without extension) | all in case all the gff 
						have the same attribute pattern
					(ii) child_attribute: attribute field of the mRNA or transcript feature in the file like ID
					(iii) child_parent_linker: attribute field of the mRNA or transcript, CDS, exon features that link them with their 
					respective parent feature like Parent - Note: base assumption by the tool is that all child levels
					have the same child-parent linker attribute fields. For eg., if Parent is the child-parent linker in the mRNA feature line,
					then Parent will be the child-parent linker for all other child-level feature lines in the GFF
					(iv) parent_attribute: attribute field of the gene feature like ID 
					Default if no config is given (or no matching row): ID, Parent, ID

					--process_pseudos <yes|no>
					yes (default): CDS of pseudogenes are extracted like all others
					no: CDS of pseudogenes are skipped

					--transsplice <warn|skip|process>
					Controls transcripts that look trans-spliced (CDS segments on different strands or carrying part= attributes)
					warn (default): print a warning and extract using the standard single-strand rules (sequence might be wrong)
					skip: print a warning and leave these transcripts out
					process: print a warning and extract them, ordering segments by their part= numbers and handling the strand of
							 each segment separately; transcripts without part= numbers on all segments are skipped

					--contig_boundary_violation_process <yes|no>
					Controls CDS segments whose GFF coordinates run past the end of the contig. NCBI writes genes that cross the
					origin of circular genomes (mitochondria, plastids) this way, e.g. a CDS segment at 394625-395019 on a 367808 bp contig
					yes: such segments wrap around to the start of the contig (circular contig assumed) and a warning is printed for each
					no (default): coordinates are used as given, so such segments come out empty or cut short without a warning;
								  check the summary for sequences whose length is not a multiple of 3

			Transcripts with unstranded features (strand other than + or -) are always skipped with a warning.
			CDS lines with several comma-separated parents are assigned to each of the parents.

			Output files in the output directory:
					<name>.cds.fasta - CDS sequences, one per transcript
					<name>.cds_summary.txt - counts of sequences, length not a multiple of 3, missing ATG start,
											 missing terminal stop, internal stops, non-ACGT characters, duplicate IDs
					<name>.cds_problems.tsv - IDs and lengths of the flagged sequences with the issues found
					cds_summary_all_species.tsv - one summary row per species
			"""

### --- begin imports --- ###
import os, sys, glob, re, subprocess, gzip
from operator import itemgetter


### --- end of imports --- ###

#Loading gene IDs from FASTA file
def load_sequences( multiple_fasta_file ):
	"""! @brief load candidate gene IDs from file """
	if multiple_fasta_file[-2:].lower() != 'gz':#dealing with uncompressed fasta file
		sequences = {}
		with open(multiple_fasta_file) as f:
			header = f.readline()[1:].strip()
			if " " in header:
				header = header.split(' ')[0]
				if "\t" in header:
					header = header.split('\t')[0]
			elif "\t" in header:
				header = header.split('\t')[0]
			seq = []
			line = f.readline()
			while line:
				if line[0] == '>':
					sequences.update({header: "".join(seq).upper()})
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
						if "\t" in header:
							header = header.split('\t')[0]
					elif "\t" in header:
						header = header.split('\t')[0]
					seq = []
				else:
					seq.append(line.strip())
				line = f.readline()
			sequences.update({header: "".join(seq).upper()})
		return sequences

	else:#dealing with compressed FASTA file(s)
		sequences = {}
		with gzip.open (multiple_fasta_file, "rt") as f:
			header = f.readline()[1:].strip()
			if " " in header:
				header = header.split(' ')[0]
				if "\t" in header:
					header = header.split('\t')[0]
			elif "\t" in header:
				header = header.split('\t')[0]
			seq = []
			line = f.readline()
			while line:
				if line[0] == '>':
					sequences.update({header: "".join(seq).upper()})
					header = line.strip()[1:]
					if " " in header:
						header = header.split(' ')[0]
						if "\t" in header:
							header = header.split('\t')[0]
					elif "\t" in header:
						header = header.split('\t')[0]
					seq = []
				else:
					seq.append(line.strip())
				line = f.readline()
			sequences.update({header: "".join(seq).upper()})
		return sequences

#function for loading transcript information from GFF3 file
def load_transcript_information_from_gff3( gff3_input_file,process_pseudos,child_attribute,child_parent_linker, parent_attribute, transsplice):
	"""! @brief load all transcript information from gff3 file """
	# --- load all data from file --- #
	message = []
	gff_pseudos = set()

	if gff3_input_file[-2:].lower() != 'gz':  # dealing with uncompressed gff file
		information = []
		mrna_dict = {}  # To store mRNA ID -> details mapping
		exon_dict = {}  # To store exon parent -> list of exons mapping
		#code to check for pseudogenes
		with open(gff3_input_file, "r") as f:
			gff_lines = f.readlines()
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines
					  if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines
					  if not line.startswith('#') and len(line.split('\t')) >= 3)
			# checking if cds feature is present in the GFF file
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines
						  if not line.startswith('#') and len(line.split('\t')) >= 3)
		if not has_cds:#early return to avoid writing output for cases where GFF files miss CDS features
			print(f"WARNING: {gff3_input_file}: no CDS features found; skipping "
				  f"(exon-based output would include UTRs, not CDS).")
			return [], message
		if process_pseudos == 'no':
			if has_mrna:  # checking for pseudogenes when mrna feature is present
				gff_pseudos_genes = set()
				with open(gff3_input_file, "r") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE':
								for attr in parts[8].split(';'):
									if attr.startswith((parent_attribute)+'='):
										gene_id = attr[(len(parent_attribute)+1):]
										gff_pseudos_genes.add(gene_id)
										break
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'MRNA':
								mrna_id = None
								parent = None
								# Parse all attributes first
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										mrna_id = attr[(len(child_attribute)+1):]
									elif attr.startswith(str(child_parent_linker)+'='):
										parent = attr[(len(child_parent_linker)+1):]
								# Check if this mRNA belongs to a pseudogene
								if parent and parent in gff_pseudos_genes:
									if mrna_id:
										gff_pseudos.add(mrna_id)
								else:
									if 'PSEUDO=TRUE' in parts[8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[8].upper():
										if mrna_id:
											gff_pseudos.add(mrna_id)
			elif has_transcript:  # checking for pseudogenes when transcript feature is present
				gff_pseudos_genes = set()
				with open(gff3_input_file, "r") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE':
								for attr in parts[8].split(';'):
									if attr.startswith(str(parent_attribute)+'='):
										gene_id = attr[(len(parent_attribute)+1):]
										gff_pseudos_genes.add(gene_id)
										break

					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'TRANSCRIPT':
								mrna_id = None
								parent = None
								# Parse all attributes first
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										mrna_id = attr[(len(child_attribute)+1):]
									elif attr.startswith(str(child_parent_linker)+'='):
										parent = attr[(len(child_parent_linker)+1):]
								# Check if this mRNA belongs to a pseudogene
								if parent and parent in gff_pseudos_genes:
									if mrna_id:
										gff_pseudos.add(mrna_id)
								else:
									if 'PSEUDO=TRUE' in parts[8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[8].upper():
										if mrna_id:
											gff_pseudos.add(mrna_id)
			else:  # checking for pseudogenes when no mrna, transcript feature is present and only cds feature is present
				with open(gff3_input_file, "r") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and (parts[2].upper() == 'PSEUDOGENE' or 'PSEUDO=TRUE' in parts[8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[8].upper()):
								for attr in parts[8].split(';'):
									if attr.startswith(str(parent_attribute)+'='):
										gene_id = attr[(len(parent_attribute)+1):]
										gff_pseudos.add(gene_id)
										break

		with open(gff3_input_file, "r") as f:
			line = f.readline()
			while line:
				if line[0] != '#':
					parts = line.strip().split('\t')
					if len(parts) > 5:
						if has_cds:
							if parts[2].upper() == 'CDS':
								if len(parts[-1]) > len(str(child_parent_linker)+'='):
									if ";" in parts[-1]:
										parent = None  # Changed from False to None
										part_value = None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if (str(child_parent_linker)+'=') in subp:
												parent = subp.replace((str(child_parent_linker)+'='), "")
											if subp.startswith('part='):#to handle trans-splicing events
												part_value = int(subp.replace('part=', '')) if subp.replace('part=','').isdigit() else None

										# Check if parent is pseudogene AFTER finding parent
										if parent:
											for each_parent in [p.strip() for p in parent.split(',')]:#to handle CDS with multiple parents
												if each_parent in gff_pseudos:
													continue # Skip this CDS

												information.append({
													'chr': parts[0],
													'start': int(parts[3]),
													'end': int(parts[4]),
													'orientation': parts[6],
													'parent': each_parent,
													'phase': int(parts[7]) if parts[7].isdigit() else 0,
													'part': part_value
												})
										else:
											message.append("no parent detected - " + line)
									else:
										parent = None
										if (str(child_parent_linker)+'=') in parts[-1]:
											parent = str(parts[-1]).replace((str(child_parent_linker)+'='), "")

										if parent:
											for each_parent in [p.strip() for p in parent.split(',')]:#to handle CDS with multiple parents
												if each_parent in gff_pseudos:
													continue #skip CDS

												information.append({
													'chr': parts[0],
													'start': int(parts[3]),
													'end': int(parts[4]),
													'orientation': parts[6],
													'parent': each_parent,
													'phase': int(parts[7]) if parts[7].isdigit() else 0,
													'part': None
												})
										else:
											message.append("only one field - " + line)
						else:
							# Handle MRNA/TRANSCRIPT and EXON features
							if parts[2].upper() == 'MRNA' or parts[2].upper() == 'TRANSCRIPT':
								if len(parts[-1]) > len((str(child_attribute)+'=')):
									if ";" in parts[-1]:
										mrna_id = None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if (str(child_attribute)+'=') in subp:
												mrna_id = subp.replace((str(child_attribute)+'='), "")

										if mrna_id and mrna_id not in gff_pseudos:
											mrna_dict[mrna_id] = {
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
											}

							elif parts[2].upper() == 'EXON':
								if len(parts[-1]) > len(str(child_parent_linker)+'='):
									if ";" in parts[-1]:
										parent = None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if (str(child_parent_linker)+'=') in subp:
												parent = subp.replace((str(child_parent_linker)+'='), "")

										if parent and parent not in gff_pseudos:
											if parent not in exon_dict:
												exon_dict[parent] = []
											exon_dict[parent].append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent,
											})
										else:
											message.append("no parent detected - " + line)
									else:
										parent = None
										if (str(child_parent_linker)+'=') in parts[-1]:
											parent = str(parts[-1]).replace((str(child_parent_linker)+'='), "")

										if parent and parent not in gff_pseudos:
											if parent not in exon_dict:
												exon_dict[parent] = []
											exon_dict[parent].append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("only one field - " + line)
				line = f.readline()

		if has_cds:
			sorted_data = {}
			for each in information:
				try:
					sorted_data[each['parent']].append(each)
				except KeyError:
					sorted_data.update({each['parent']: [each]})

			final_data = []
			for key in sorted_data.keys():
				segments = sorted_data[key]
				orientations = set(s['orientation'] for s in segments)
				has_part = any(s.get('part') is not None for s in segments)

				if orientations - {'+','-'}:#skipping transcripts that do not have + or - orientation
					print(f"WARNING: transcript '{key}' has unstranded (.) features; skipped.")  # NEW
					continue

				if len(orientations) > 1 or has_part:
					print(f"WARNING: transcript '{key}' shows signs of trans-splicing "
								   f"(mixed strands: {len(orientations) > 1}, part= present: {has_part}) — "
								   f"--transsplice mode: {transsplice}.")
					if transsplice == 'skip':
						continue
					elif transsplice == 'process':
						if all(s['part'] is not None for s in segments):
							final_data.append(sorted(segments, key=lambda s: (s['part'], s['start'] if s['orientation'] == '+' else -s['start'])))
						else:
							print(f"WARNING: transcript '{key}' lacks part= numbers on some segments; segment order cannot be determined, skipped.")
						continue
				if segments[0]['orientation'] == '+':
					final_data.append(sorted(segments, key=itemgetter('start')))
				else:
					final_data.append(sorted(segments, key=itemgetter('start'))[::-1])
		else:
			# Process exon/transcript based annotation
			final_data = []
			for mrna_id, mrna_info in mrna_dict.items():
				if mrna_id in exon_dict:
					exons = exon_dict[mrna_id]

					# Sort the exons based on orientation
					if mrna_info['orientation'] == '+':
						sorted_exons = sorted(exons, key=itemgetter('start'))
					else:
						sorted_exons = sorted(exons, key=itemgetter('start'), reverse=True)
					final_data.append(sorted_exons)
		return final_data, message

	else:#dealing with compressed GFF3 file(s)
		information = []
		mrna_dict = {}  # To store mRNA ID -> details mapping
		exon_dict = {}  # To store exon parent -> list of exons mapping


		with gzip.open(gff3_input_file, "rt") as f:
			gff_lines = f.readlines()
			has_mrna = any(line.split('\t')[2].upper() == 'MRNA' for line in gff_lines
				   if not line.startswith('#') and len(line.split('\t')) >= 3)
			has_transcript = any(line.split('\t')[2].upper() == 'TRANSCRIPT' for line in gff_lines
						 if not line.startswith('#') and len(line.split('\t')) >= 3)
			# checking if cds feature is present in the GFF file
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines
				  if not line.startswith('#') and len(line.split('\t')) >= 3)
		if not has_cds:#early return to avoid writing output for cases where GFF files miss CDS features
			print(f"WARNING: {gff3_input_file}: no CDS features found; skipping "
				  f"(exon-based output would include UTRs, not CDS).")
			return [], message

		"""
		with gzip.open(gff3_input_file, "rt") as f:
			gff_lines = f.readlines()
			# collecting pseudogenes to skip cds of pseudogenes
			for line in gff_lines:
				if line[0] != '#':
					parts = line.strip().split('\t')
					if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE':
						for attr in parts[8].split(';'):
							if attr.startswith(str(child_attribute)+'='):
								gene_id = attr[len(str(child_attribute))+1:]
								gff_pseudos.add(gene_id)
								break

			# checking if cds feature is present in the GFF file
			has_cds = any(line.split('\t')[2].upper() == 'CDS' for line in gff_lines
						  if not line.startswith('#') and len(line.split('\t')) >= 3)
		"""

		if process_pseudos == 'no':
			if has_mrna:  # checking for pseudogenes when mrna feature is present
				gff_pseudos_genes = set()
				with gzip.open(gff3_input_file, "rt") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE':
								for attr in parts[8].split(';'):
									if attr.startswith(str(parent_attribute)+'='):
										gene_id = attr[len(str(parent_attribute))+1:]
										gff_pseudos_genes.add(gene_id)
										break
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'MRNA':
								mrna_id = None
								parent = None
								# Parse all attributes first
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										mrna_id = attr[len(str(child_attribute))+1:]
									elif attr.startswith(str(child_parent_linker)+'='):
										parent = attr[len(str(child_parent_linker))+1:]
								# Check if this mRNA belongs to a pseudogene
								if parent and parent in gff_pseudos_genes:
									if mrna_id:
										gff_pseudos.add(mrna_id)
								else:
									if 'PSEUDO=TRUE' in parts[8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[8].upper():
										if mrna_id:
											gff_pseudos.add(mrna_id)
			elif has_transcript:  # checking for pseudogenes when transcript feature is present
				gff_pseudos_genes = set()
				with gzip.open(gff3_input_file, "rt") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'PSEUDOGENE':
								for attr in parts[8].split(';'):
									if attr.startswith(str(parent_attribute)+'='):
										gene_id = attr[len(str(parent_attribute))+1:]
										gff_pseudos_genes.add(gene_id)
										break

					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and parts[2].upper() == 'TRANSCRIPT':
								mrna_id = None
								parent = None
								# Parse all attributes first
								for attr in parts[8].split(';'):
									if attr.startswith(str(child_attribute)+'='):
										mrna_id = attr[len(str(child_attribute))+1:]
									elif attr.startswith(str(child_parent_linker)+'='):
										parent = attr[len(str(child_parent_linker))+1:]
								# Check if this mRNA belongs to a pseudogene
								if parent and parent in gff_pseudos_genes:
									if mrna_id:
										gff_pseudos.add(mrna_id)
								else:
									if 'PSEUDO=TRUE' in parts[8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[8].upper():
										if mrna_id:
											gff_pseudos.add(mrna_id)
			else:  # checking for pseudogenes when no mrna, transcript feature is present and only cds feature is present
				with gzip.open(gff3_input_file, "rt") as f:
					gff_lines = f.readlines()
					# collecting pseudogenes to skip cds of pseudogenes
					for line in gff_lines:
						if line[0] != '#':
							parts = line.strip().split('\t')
							if len(parts) >= 9 and (parts[2].upper() == 'PSEUDOGENE' or 'PSEUDO=TRUE' in parts[8].upper() or 'GENE_BIOTYPE=PSEUDOGENE' in parts[8].upper()):
								for attr in parts[8].split(';'):
									if attr.startswith(str(parent_attribute)+'='):
										gene_id = attr[len(str(parent_attribute))+1:]
										gff_pseudos.add(gene_id)
										break

		with gzip.open(gff3_input_file, "rt") as f:
			line = f.readline()
			while line:
				if line[0] != '#':
					parts = line.strip().split('\t')
					if len(parts) > 5:
						if has_cds:
							if parts[2].upper() == 'CDS':
								if len(parts[-1]) > len(str(child_parent_linker)+'='):
									if ";" in parts[-1]:
										parent = None  # Changed from False to None
										part_value = None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if str(child_parent_linker)+'=' in subp:
												parent = subp.replace(str(child_parent_linker)+'=', "")
											if subp.startswith('part='):#to handle trans-splicing events
												part_value = int(subp.replace('part=', '')) if subp.replace('part=','').isdigit() else None

										# Check if parent is pseudogene AFTER finding parent
										if parent:
											for each_parent in [p.strip() for p in parent.split(',')]:#to handle CDS with multiple parents
												if each_parent in gff_pseudos:
													continue

												information.append({
													'chr': parts[0],
													'start': int(parts[3]),
													'end': int(parts[4]),
													'orientation': parts[6],
													'parent': each_parent,
													'phase': int(parts[7]) if parts[7].isdigit() else 0,
													'part': part_value
												})
										else:
											message.append("no parent detected - " + line)
									else:
										parent = None
										if str(child_parent_linker)+'=' in parts[-1]:
											parent = str(parts[-1]).replace(str(child_parent_linker)+'=', "")

										if parent:
											for each_parent in [p.strip() for p in parent.split(',')]:#to handle CDS with multiple parents
												if each_parent in gff_pseudos:
													continue#skip CDS

												information.append({
													'chr': parts[0],
													'start': int(parts[3]),
													'end': int(parts[4]),
													'orientation': parts[6],
													'parent': each_parent,
													'phase': int(parts[7]) if parts[7].isdigit() else 0,
													'part': None
												})
										else:
											message.append("only one field - " + line)
						else:
							# Handle MRNA/TRANSCRIPT and EXON features
							if parts[2].upper() == 'MRNA' or parts[2].upper() == 'TRANSCRIPT':
								if len(parts[-1]) > len(str(child_attribute)+'='):
									if ";" in parts[-1]:
										mrna_id = None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if str(child_attribute)+'=' in subp:
												mrna_id = subp.replace(str(child_attribute)+'=', "")

										if mrna_id and mrna_id not in gff_pseudos:
											mrna_dict[mrna_id] = {
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
											}

							elif parts[2].upper() == 'EXON':
								if len(parts[-1]) > len(str(child_parent_linker)+'='):
									if ";" in parts[-1]:
										parent = None
										subparts = parts[-1].split(';')
										for subp in subparts:
											if str(child_parent_linker)+'=' in subp:
												parent = subp.replace(str(child_parent_linker)+'=', "")

										if parent and parent not in gff_pseudos:
											if parent not in exon_dict:
												exon_dict[parent] = []
											exon_dict[parent].append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("no parent detected - " + line)
									else:
										parent = None
										if str(child_parent_linker)+'=' in parts[-1]:
											parent = str(parts[-1]).replace(str(child_parent_linker)+'=', "")

										if parent and parent not in gff_pseudos:
											if parent not in exon_dict:
												exon_dict[parent] = []
											exon_dict[parent].append({
												'chr': parts[0],
												'start': int(parts[3]),
												'end': int(parts[4]),
												'orientation': parts[6],
												'parent': parent
											})
										else:
											message.append("only one field - " + line)
				line = f.readline()

		if has_cds:
			# --- sort data by parent --- #
			sorted_data = {}
			for each in information:
				try:
					sorted_data[each['parent']].append(each)
				except KeyError:
					sorted_data.update({each['parent']: [each]})

			final_data = []
			for key in sorted_data.keys():
				segments = sorted_data[key]
				orientations = set(s['orientation'] for s in segments)
				has_part = any(s.get('part') is not None for s in segments)

				if orientations - {'+','-'}:#skipping transcripts that do not have + or - orientation
					print(f"WARNING: transcript '{key}' has unstranded (.) features; skipped.")  # NEW
					continue

				if len(orientations) > 1 or has_part:
					print(f"WARNING: transcript '{key}' shows signs of trans-splicing "
								   f"(mixed strands: {len(orientations) > 1}, part= present: {has_part}) — "
								   f"--transsplice mode: {transsplice}.")
					if transsplice == 'skip':
						continue
					elif transsplice == 'process':
						if all(s['part'] is not None for s in segments):
							final_data.append(sorted(segments, key=lambda s: (s['part'], s['start'] if s['orientation'] == '+' else -s['start'])))
						else:
							print(f"WARNING: transcript '{key}' lacks part= numbers on some segments; segment order cannot be determined, skipped.")
						continue

				if segments[0]['orientation'] == '+':
					final_data.append(sorted(segments, key=itemgetter('start')))
				else:
					final_data.append(sorted(segments, key=itemgetter('start'))[::-1])
		else:
			# Process exon/transcript based annotation
			final_data = []
			for mrna_id, mrna_info in mrna_dict.items():
				if mrna_id in exon_dict:
					exons = exon_dict[mrna_id]
					# Sort the exons based on orientation
					if mrna_info['orientation'] == '+':
						sorted_exons = sorted(exons, key=itemgetter('start'))
					else:
						sorted_exons = sorted(exons, key=itemgetter('start'), reverse=True)
					final_data.append(sorted_exons)
		return final_data, message

#function of constructing the CDS FASTA file
def construct_CDS_file( transcript_info, CDS_file, assembly, child_parent_linker, transsplice, process_violating_segments):
	"""! @brief construct file with all sequences for translation """
	with open(CDS_file, "w") as out:
		for transcript in transcript_info:
			if any(p['chr'] not in assembly for p in transcript):
				continue
			seq = []
			revcomp_status = False
			if transcript[0]['orientation'] == '-':
				revcomp_status = True
			for part in transcript:
				if transsplice == 'process':#to handle transspliced segments
					revcomp_status = part['orientation'] == '-'
				if revcomp_status:
					if process_violating_segments=='yes':
						seq.append(revcomp(get_segment(assembly[part['chr']],part['start'],part['end'])))
					else:
						seq.append(revcomp(assembly[part['chr']][part['start'] - 1:part['end']]))
				else:
					if process_violating_segments=='yes':
						seq.append(get_segment(assembly[part['chr']],part['start'],part['end']))
					else:
						seq.append(assembly[part['chr']][part['start'] - 1:part['end']])
				# Get parent ID, handling both formats of CDS feature or mRNA/ exon features' presence
				parent_id = transcript[0]['parent']
				if str(child_parent_linker) + "=" in parent_id:
					parent_id = parent_id.replace(str(child_parent_linker) + "=", "")
			first_phase = transcript[0].get('phase', 0)#to account for different frames (phase values) for slicing the sequences
			if first_phase:
				seq[0] = seq[0][first_phase:]
			out.write( '>' + str(parent_id) + '\n' + "".join( seq ) + '\n' )

#function for constructing the reverse complement
def revcomp( seq ):
	"""! @brief constructs revcomp """
	new_seq = []
	dictionary = { 'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N',
		'R':'Y', 'Y':'R', 'S':'S', 'W':'W',
		'K':'M', 'M':'K', 'B':'V', 'V':'B', 'D':'H', 'H':'D',
		'a':'t', 't':'a', 'c':'g', 'g':'c', 'n':'n',
		'r':'y', 'y':'r', 's':'s', 'w':'w',
		'k':'m', 'm':'k', 'b':'v', 'v':'b', 'd':'h', 'h':'d' }#modified dictionary to handle all IUPAC codes
	for nt in seq:
		try:
			new_seq.append( dictionary[ nt ] )
		except KeyError:
			print("WARNING: unexpected character in sequence during revcomp: '" + nt + "'")
			new_seq.append( "N")
	return ''.join( new_seq[::-1] )

#function to return basename without extension
def strip_extension(file_name, extensions):
	"""returns the base name without the extension, or None if no extension matches"""
	lower_name = file_name.lower()
	for ext in extensions:
		if lower_name.endswith(ext):
			return lower_name[:-len(ext)]
	return None

#function to summarize the output CDS file and generate a summary report of a CDS FASTA file
def summarize_cds_file( cds_file, report_file, problem_file ):
	"""! @brief writes a summary report and a list of problematic sequences for one CDS FASTA file """
	stop_codons = ['TAA', 'TAG', 'TGA']

	#reads all records as a list so that duplicate IDs are not lost
	records = []
	with open(cds_file) as f:
		header = None
		seq = []
		for line in f:
			line = line.strip()
			if line.startswith('>'):
				if header is not None:
					records.append((header, "".join(seq)))
				header = line[1:]
				seq = []
			else:
				seq.append(line)
		if header is not None:
			records.append((header, "".join(seq)))

	total = len(records)
	not_div_by_3 = 0
	no_start = 0
	no_stop = 0
	internal_stop = 0
	non_acgt = 0
	duplicate_ids = 0
	clean = 0
	lengths = []
	seen_ids = set()

	with open(problem_file, "w") as prob:
		prob.write("id\tlength\tissues\n")
		for header, seq in records:
			issues = []
			length = len(seq)
			lengths.append(length)

			if length % 3 != 0:
				not_div_by_3 += 1
				issues.append("length_not_multiple_of_3")
			if seq[:3] != 'ATG':
				no_start += 1
				issues.append("no_ATG_start")
			if seq[-3:] not in stop_codons:
				no_stop += 1
				issues.append("no_terminal_stop")

			# in-frame codons (incomplete last codon ignored); last codon excluded because a terminal stop is expected
			body = seq[:-3] if (length % 3 == 0 and seq[-3:] in stop_codons) else seq
			codons = [body[i:i + 3] for i in range(0, len(body) - 2, 3)]
			if any(codon in stop_codons for codon in codons):
				internal_stop += 1
				issues.append("internal_stop")

			if any(nt not in 'ACGT' for nt in seq):
				non_acgt += 1
				issues.append("non_ACGT_characters")

			if header in seen_ids:
				duplicate_ids += 1
				issues.append("duplicate_id")
			seen_ids.add(header)

			if issues:
				prob.write(header + "\t" + str(length) + "\t" + ",".join(issues) + "\n")
			else:
				clean += 1

	def pct(n):
		if total == 0:
			return "NA"
		return f"{100.0 * n / total:.2f}%"

	lines = [
		f"CDS file: {cds_file}",
		f"Total sequences written: {total}",
		f"Clean (multiple of 3, ATG start, terminal stop, no internal stop): {clean} ({pct(clean)})",
		f"Length not a multiple of 3: {not_div_by_3} ({pct(not_div_by_3)})",
		f"No ATG start: {no_start} ({pct(no_start)})",
		f"No terminal stop codon: {no_stop} ({pct(no_stop)})",
		f"Internal stop codon(s): {internal_stop} ({pct(internal_stop)})",
		f"Non-ACGT characters (e.g. N): {non_acgt} ({pct(non_acgt)})",
		f"Duplicate sequence IDs: {duplicate_ids}",
	]
	if lengths:
		sorted_lengths = sorted(lengths)
		lines.append(f"Length min / median / mean / max: {sorted_lengths[0]} / {sorted_lengths[len(sorted_lengths)//2]} / {sum(lengths)/len(lengths):.1f} / {sorted_lengths[-1]}")
	report = "\n".join(lines) + "\n"
	with open(report_file, "w") as out:
		out.write(report)
	print(report)
	return [total, clean, not_div_by_3, no_start, no_stop, internal_stop, non_acgt, duplicate_ids]

#function to extract a segment, handling segments that run past the end of a circular contig
def get_segment( contig_seq, start, end ):
	"""! @brief returns the sequence of a segment; coordinates beyond the contig length wrap around to the start """
	contig_length = len(contig_seq)
	if end <= contig_length:
		return contig_seq[start - 1:end]
	print(f"WARNING: segment {start}-{end} extends beyond contig length {contig_length}; wrapping around (circular contig assumed)")
	if start > contig_length:#whole segment lies beyond the end
		return contig_seq[start - contig_length - 1:end - contig_length]
	return contig_seq[start - 1:] + contig_seq[:end - contig_length]#segment crosses the end

def main (arguments):
	print(f'You are using get_cds_from_assembly.py version {__version__}.')
	assembly_exts = ['.fasta.gz', '.fna.gz', '.fa.gz', '.fasta', '.fna', '.fa']
	gff_exts = ['.gff3.gz', '.gff.gz', '.gff3', '.gff']
	assembly_arg=arguments[arguments.index('--assembly')+1]#full path to file or folder of assembly files
	gff_arg=arguments[arguments.index('--gff')+1]#full path to file or folder of gff files
	outdir = arguments[arguments.index('--out') + 1]
	if outdir[-1] != "/":
		outdir += "/"
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	if '--process_pseudos' in arguments:#yes/no
		process_pseudos = arguments[arguments.index('--process_pseudos') + 1]
	else:
		process_pseudos = 'yes'
	if '--transsplice' in arguments:#warn throws a warning and processes the segments with standard orientation rules (might be wrong), skip skips the transspliced segments entirely, process helps process the transspliced segments according to the parts numbering in case they are present
		transsplice = arguments[arguments.index('--transsplice') + 1].lower()
		if transsplice not in ('warn', 'skip', 'process'):
			sys.exit("--transsplice must be 'warn', 'skip' or 'process'")
	else:
		transsplice = 'warn'

	if '--contig_boundary_violation_process' in arguments:# yes or no to process segments with GFF boundaries beyond the contig length, usual for circular contigs like mitochondria
		process_violating_segments = arguments[arguments.index('--contig_boundary_violation_process') + 1]
	else:
		process_violating_segments = 'no'
	gff_config_parameters = {}
	if '--gff_config' in arguments:
		gff_config_file = arguments[arguments.index('--gff_config')+1]
		with open(gff_config_file,'r')as f:
			for line in f:
				parts=line.strip().split()
				if len(parts)==4:
					gff_config_parameters[(parts[0].lower())] = {'child_attribute': parts[1],'child_parent_linker': parts[2],'parent_attribute': parts[3]}
	else:
		gff_config_parameters['default'] = {
			'child_attribute': 'ID',
			'child_parent_linker': 'Parent',
			'parent_attribute': 'ID'
		}

	assembly_files=[]
	gff_files=[]
	if os.path.isfile(assembly_arg):
		assembly_files=[assembly_arg]
	elif os.path.isdir(assembly_arg):
		for f in sorted(glob.glob(os.path.join(assembly_arg, "*"))):
			if strip_extension(os.path.basename(f), assembly_exts) is not None:
				assembly_files.append(f)
	if os.path.isfile(gff_arg):
		gff_files=[gff_arg]
	elif os.path.isdir(gff_arg):
		for f in sorted(glob.glob(os.path.join(gff_arg, "*"))):
			if strip_extension(os.path.basename(f), gff_exts) is not None:
				gff_files.append(f)
	#code block to check for missing pair of gff or assembly file(s)
	a_names = [strip_extension(os.path.basename(a), assembly_exts) for a in assembly_files]
	g_names = [strip_extension(os.path.basename(g), gff_exts) for g in gff_files]
	for n in a_names:
		if n not in g_names:
			print(f"WARNING: assembly '{n}' has no matching GFF")
	for n in g_names:
		if n not in a_names:
			print(f"WARNING: GFF '{n}' has no matching assembly")
	all_stats = []
	for assembly in assembly_files:
		assembly_name = strip_extension(os.path.basename(assembly), assembly_exts)
		if assembly_name is None:
			print(f"WARNING: {assembly} does not have a recognised assembly extension, skipped")
			continue
		for gff in gff_files:
			gff_name = strip_extension(os.path.basename(gff), gff_exts)
			if gff_name is None:
				print(f"WARNING: {gff} does not have a recognised GFF extension, skipped")
				continue
			if assembly_name == gff_name:
				cds_out = os.path.join(outdir,f"{gff_name}.cds.fasta")
				print(f"gff_name is {gff_name}")
				if gff_name in gff_config_parameters:
					child_attribute = gff_config_parameters[gff_name]['child_attribute']
					child_parent_linker = gff_config_parameters[gff_name]['child_parent_linker']
					parent_attribute = gff_config_parameters[gff_name]['parent_attribute']
				elif 'all' in gff_config_parameters:
					child_attribute = gff_config_parameters['all']['child_attribute']
					child_parent_linker = gff_config_parameters['all']['child_parent_linker']
					parent_attribute = gff_config_parameters['all']['parent_attribute']
				else:
					child_attribute = gff_config_parameters['default']['child_attribute']
					child_parent_linker = gff_config_parameters['default']['child_parent_linker']
					parent_attribute = gff_config_parameters['default']['parent_attribute']
				seqs = load_sequences(str(assembly))
				transcript_information, message = load_transcript_information_from_gff3(str(gff), process_pseudos,child_attribute,child_parent_linker, parent_attribute, transsplice)
				if not transcript_information:#check to skip the rest of the code block in case an empty transcript_information is returned from the load_transcript_information function
					print(f"WARNING: {gff_name}: no transcripts extracted, no output written")
					continue
				#warning in case GFF contig names do not match with the FASTA header contig names
				gff_ids = {seg['chr'] for transcript in transcript_information for seg in transcript}
				missing_ids = gff_ids - set(seqs)
				if gff_ids - set(seqs):
					print(f"WARNING: {gff_name}: {len(missing_ids)} of {len(gff_ids)} GFF seqids not in FASTA, e.g. {sorted(missing_ids)[:5]}")
					if len(missing_ids) == len(gff_ids):
						print(f"ERROR: {gff_name}: no GFF seqid matches the FASTA; skipping this species")
						continue
				construct_CDS_file(transcript_information, cds_out, seqs, child_parent_linker, transsplice,process_violating_segments)
				stats = summarize_cds_file(cds_out, os.path.join(outdir, f"{gff_name}.cds_summary.txt"),os.path.join(outdir, f"{gff_name}.cds_problems.tsv"))
				all_stats.append([gff_name] + stats)
	if all_stats:
		with open(os.path.join(outdir, "cds_summary_all_species.tsv"), "w") as out:
			out.write("species\ttotal\tclean\tlength_not_multiple_of_3\tno_ATG_start\tno_terminal_stop\tinternal_stop\tnon_ACGT\tduplicate_ids\n")
			for row in all_stats:
				out.write("\t".join(str(x) for x in row) + "\n")

if '--assembly' in sys.argv and '--gff' in sys.argv and '--out' in sys.argv:
	main(sys.argv)
else:
	sys.exit(__usage__)