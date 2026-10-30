### V0.1 ###
### Shakunthala Natarajan ###
"""
feature requests and bug reports to s64snata@uni-bonn.de
"""

__version__='0.1'
__usage__= """
			python3 compare_cds_seq.py
			--reference <full path to folder of reference CDS files>
			--test <full path to folder of CDS files to be tested>
			--out <full path to output directory>

			optional:
			--seqkit <full path to the seqkit executable> (default: seqkit, i.e. found via PATH)

			NOTE_1: CDS files should have extension of .cds.fa, .cds.fna or .cds.fasta (optionally followed by .gz)
			NOTE_2: A reference file and a test file are paired when their names are identical after removing
					the extension (case ignored), e.g. species1.cds.fasta and Species1.cds.fa.gz
			NOTE_3: The comparison is done with respect to the reference. Sequences are compared by the hash
					calculated by seqkit (seqkit fx2tab -s), sequence IDs are not used for the matching
			NOTE_4: seqkit has to be installed (conda install -c bioconda seqkit). If it is not available in PATH,
					give the location of the executable with --seqkit

			Output files in the output directory:
					<name>_differences.tsv - sequences only found in the reference or only found in the test file
											 (only written if there are differences)
					comparison_summary.tsv - one summary row per name
			"""

### --- begin imports --- ###
import os, sys, shutil, subprocess


### --- end of imports --- ###

#longer extensions first, so that .cds.fa.gz is removed completely and not only .gz
EXTENSIONS = ['.cds.fasta.gz', '.cds.fna.gz', '.cds.fa.gz', '.cds.fasta', '.cds.fna', '.cds.fa']

#function to return basename without extension
def strip_extension( file_name, extensions ):
	"""! @brief returns the lowercase base name without the extension, or None if no extension matches """
	lower_name = file_name.lower()
	for ext in extensions:
		if lower_name.endswith(ext):
			return lower_name[:-len(ext)]
	return None

#function to collect the CDS files of a folder
def find_cds_files( folder ):
	"""! @brief returns a dictionary: base name (lowercase) -> full path of the CDS file """
	files = {}
	for file_name in sorted(os.listdir(folder)):
		base = strip_extension(file_name, EXTENSIONS)
		if base is None:
			continue
		if base in files:
			print(f"WARNING: {folder}: more than one file for '{base}', using {os.path.basename(files[base])}")
		else:
			files[base] = os.path.join(folder, file_name)
	return files

#function to get the sequence IDs and hashes with seqkit
def get_records( fasta_file,seqkit ):
	"""! @brief runs seqkit and returns a list of (sequence ID, sequence hash) for every sequence in the file """
	result = subprocess.run([seqkit, 'fx2tab', '-s', fasta_file], capture_output=True, text=True)
	if result.returncode != 0:
		print(f"WARNING: seqkit failed on {fasta_file}: {result.stderr.strip()}")
		return None
	records = []
	for line in result.stdout.splitlines():
		if not line.strip():
			continue
		columns = line.split('\t')
		records.append((columns[0], columns[-1]))#first column = ID, last column = hash
	return records

#function to count how often each hash occurs
def count_hashes( records ):
	"""! @brief returns a dictionary: hash -> number of occurrences """
	counts = {}
	for seq_id, seq_hash in records:
		counts[seq_hash] = counts.get(seq_hash, 0) + 1
	return counts

#function to find the sequences that are missing in the other file
def find_unmatched( records, other_records ):
	"""! @brief returns the records whose sequence cannot be found in other_records
	(a sequence that occurs twice needs to occur twice in the other file as well) """
	available = count_hashes(other_records)
	unmatched = []
	for seq_id, seq_hash in records:
		if available.get(seq_hash, 0) > 0:
			available[seq_hash] -= 1
		else:
			unmatched.append((seq_id, seq_hash))
	return unmatched

def main (arguments):
	print(f'You are using compare_cds_seq.py version {__version__}.')
	reference_arg=arguments[arguments.index('--reference')+1]#full path to folder of reference CDS files
	test_arg=arguments[arguments.index('--test')+1]#full path to folder of CDS files to be tested
	outdir = arguments[arguments.index('--out') + 1]#full path to output directory
	if outdir[-1] != "/":
		outdir += "/"
	if '--seqkit' in arguments:
		seqkit = arguments[arguments.index('--seqkit') + 1]
	else:
		seqkit = 'seqkit'
	if not os.path.isdir(reference_arg):
		sys.exit(f"{reference_arg} is not a folder")
	if not os.path.isdir(test_arg):
		sys.exit(f"{test_arg} is not a folder")
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	reference_files = find_cds_files(reference_arg)
	test_files = find_cds_files(test_arg)
	if not reference_files:
		sys.exit(f"No CDS files with an allowed extension found in {reference_arg}")

	#files that have no partner in the other folder
	for base in test_files:
		if base not in reference_files:
			print(f"WARNING: '{base}' is in the test folder but has no reference file")

	summary_rows = []
	for base, reference_file in reference_files.items():
		if base not in test_files:
			print(f"WARNING: no test file found for reference '{base}'")
			continue
		reference_records = get_records(reference_file,seqkit)
		test_records = get_records(test_files[base],seqkit)
		if reference_records is None or test_records is None:
			continue

		#reference sequences missing in the test file, and test sequences that are not in the reference
		reference_only = find_unmatched(reference_records, test_records)
		test_only = find_unmatched(test_records, reference_records)
		identical = len(reference_records) - len(reference_only)

		summary_rows.append([base, len(reference_records), len(test_records), identical, len(reference_only), len(test_only)])
		print(f"{base}: reference={len(reference_records)} test={len(test_records)} identical={identical} reference_only={len(reference_only)} test_only={len(test_only)}")

		#one file per name listing the sequences that differ (only written if there are differences)
		if reference_only or test_only:
			with open(os.path.join(outdir, f"{base}_differences.tsv"), "w") as out:
				out.write("status\tid\tsequence_hash\n")
				for seq_id, seq_hash in reference_only:
					out.write(f"reference_only\t{seq_id}\t{seq_hash}\n")
				for seq_id, seq_hash in test_only:
					out.write(f"test_only\t{seq_id}\t{seq_hash}\n")

	#one summary table for all names
	with open(os.path.join(outdir, "comparison_summary.tsv"), "w") as out:
		out.write("name\treference_total\ttest_total\tidentical\treference_only\ttest_only\n")
		for row in summary_rows:
			out.write("\t".join(str(x) for x in row) + "\n")

if '--reference' in sys.argv and '--test' in sys.argv and '--out' in sys.argv:
	main(sys.argv)
else:
	sys.exit(__usage__)