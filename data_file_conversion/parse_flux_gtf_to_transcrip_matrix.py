#/usr/bin/python

# The following module(s) are required for listing the files in a directory.
from os import listdir
from os.path import isfile, join

# The following module(s) are rquired for regular expressions.
import re

class ParseGTF:
	
	# Function that initializes the ParseGTF class.
	def __init__(self, pathToGTF):

		# Verbose information for the programmer.
		print("Initializing ParseGTF...\n")
		
		# Specified path to all annotated gene files.
		self.path = pathToGTF
		print("Specified path:\n{}\n".format(self.path))

		# Empty paramater for the file stream, having a global access to the file.
		self.transcriptExpressionCSV = ""

		# Counts the amount of processed files.
		self.gtf_index_count = 0

		# Saves the indexed content for each annotated gene file.
		#self.gtf_index = { file_name : [self.gtf.file.index[ element ], ... }
		self.gtf_index = {}

		# The indexed content for an annotated gene file.
		#self.gtf.file.index = { transcript_id : [ gene_id, raw_count ], ... }
		self.gtf_file_index = {}

		# A list with all the 'human readable' sample names.
		self.gtf_names = []
		self.alternate = {}		

	# Function that reads in GTF Files one by one, and closes the file stream afterward.
	def GTFReader(self):

		# Verbose information for the programmer.
		print("Reading GTF File...\n")

		# Store all files in the specified path.
		files = [f for f in listdir(self.path) if isfile(join(self.path, f))]
		
		#files = ['out.batch4_TCC-18-2_t10_150622_SN163_0652_AC7EUNACXX_L6_TGACCA.gtf', 'out.batch5_TCC-20-1_t10_160408_SN163_0708_BHMJTMBCXX_L2_TAGCTT_1.gtf']

		# Filters file by Origin based on the file name information.
		gtf_files = self.filterByFileOrigin(files)			

		# For file in filtered gtf_files.
		for file in gtf_files:

			# Add the file to self.gtf_file_index
			self.gtf_index[file] = self.gtf_file_index

			# Add the sub-dir to the file.
			file = "FluxCapacitor/"+file

			# Open file handle on variable gtf.
			gtf = open(file, "r")

			# Call GTFParser and start parsing the file.
			self.GTFParser(gtf)

			# Close the file handle, ending the stream for this partiuclar file in the loop.
			gtf.close()

	# Function that filters the files in the directory by file name.
	def filterByFileOrigin(self, files):

		# Verbose information for the programmer.
		print("Filtering by File Origin...")

		# Create an empty list to append all the 'correct' files.
		gtf_files = []

		# Compose a pattern to filter by file name, since the sample is in the file name we can easily pick out the samples we won't need.
		pattern = re.compile(r'(out.(batch\d{1}_TCC-\d{2}-\d{1}_t\d+).*)')

		# For every file in the listed list of files.
		for file in files:
			
			# Match each file against the pattern.
			m = re.match(pattern, file)

			# If regex finds a match between the pattern and the file name, thus filtering the files on file name.
			if m:
				# Increase the gtf_index_count with 1, telling us later that we have a total of n files.
				self.gtf_index_count += 1
				
				# Append the found files to the gtf_files list.
				gtf_files.append(m.group(1))

				# Append the found files to the gtf_names list.
				self.gtf_names.append(m.group(2))
		
		# return the filtered gtf files.
		return gtf_files

	# Function that parses the gtf file handle.
	def GTFParser(self, fileHandle):
		
		# Verbose information for the programmer.
		#print("Parsing GTF File...\n")

		# Read in the file handle.
		gtf_content = fileHandle.readlines()

		# For each line in the gtf file.
		for line in gtf_content:

			# If the line validates.
			if self.validateLineConsistency(line):

				# Split the line into bits, containing elements of lenght 9.
				bits = line.split("\t")

				# Split the 9th element (8) on the ";".
				bits = bits[8].split(";")

				# Save transcript id, gene id and reads.
				transcript_id = bits[0]
				gene_id = bits[2]
				reads = bits[3]
				
				# Store the transcript id as key in the gtf_file_index variable and assign gene_id and their reads as values.
				self.gtf_file_index[transcript_id] = [gene_id, reads]
				id = transcript_id + " " + gene_id
				if id not in self.alternate.keys():
					self.alternate[id] = [reads]
				else:
					self.alternate[id].append(reads)
		
	def validateLineConsistency(self, line):

		# Set a flag for the boolean.
		validity = False

		# Splits each line, generating exactly 9 elements.
		line_elements = line.split("\t")

		# Change this according to the TODO below.
		if len(line_elements) != 9:
			print("something went wrong with:\n")
			print(line_elements)
			sys.exit(1)
		else:
			validity = True

		# Returns a True or False according to the validity check.
		return validity

	def assembleTranscriptExpressionMatrix(self):
		print("Assembling Transcript Expression Matrix...\n\n")
		
		# Opens a file handle for the new Transcript Expression Matrix.
		self.createTranscriptExpressionCSV()
		
		# Store all the keys for the files.
		gtf_keys = self.gtf_index.keys()

		# Create an empty list that will contain row's.
		rows = []

		# For each key in the stored keys.
		for key in gtf_keys:

			# Save the values from the files' keys. This is another dictionary with the content of the file.
			content = self.gtf_index[key]

			#self.writeToTranscriptExpressionCSV(content)
			# Store the keys (Transcript id's) for each file.
			content_keys = content.keys()

			# For every key in the content of the file.
			for c_key in content_keys:
			
				# Save the value pair from the keys (Gene ID and read count).
				values = content[c_key]

				# Splitting the key results in the "ENSTXXXXX".
				c_k = c_key.split(" ")
				tr_id = c_k[1]
				
				# Splitting the first element of the values results in a gene id literal and the gene id. We save the third element, which is the gene_id.
				gene_id = values[0].split(" ")
				gene_id = gene_id[2]

				# Splitting the 2nd element of the values results in a reads literal and the read count. We save the third element, which is the read count.
				reads = values[1].split(" ")
				reads = reads[2]

				# Assemble row.
				row = tr_id + "\t" + gene_id + "\t" + reads + "\n" 

				# Add row to rows.
				rows.append(row)

		#self.writeToTranscriptExpressionCSV(rows)
		self.writeToTranscriptExpressionCSV()		
	def createTranscriptExpressionCSV(self):
		print("Creating Transcript Expression CSV...")

		# Creates a class wide filehandle with name gsTcell_TranscriptExpression.csv
		self.transcriptExpressionCSV = open("gsTcell_TranscriptExpression.csv","w")
		
		# Create a string with sample names.
		sample_row = "trId\tgeneId\t"
		for i in self.gtf_names:
			sample_row += i + "\t"
		sample_row = sample_row + "\n"
		
		# Write sample_row to transcript expression file. This happens only once!
		self.transcriptExpressionCSV.write(sample_row)
		

	def writeToTranscriptExpressionCSV(self):
		#print("Writing to Transcript Expression CSV...")
		
		# Store all the keys from the content of the gtf file.
		keys = self.alternate.keys()

		# Write row to file.
		#self.transcriptExpressionCSV.write()

		# For each key in the keys list.
		row = ""
		
		for key in keys:

			# Get the value pairs from the keys.
			#values = content[key]
			counts = self.alternate[key]
			#print(self.alternate)
			# Split the key, so that the transcript id string is separated from the transcript id.
			#k = key.split(" ")
			id = key.split(" ")
			
			id = id[1] + "\t" + id[4]
			
			row += id + "\t"
			# Split the gene id value, to separate the literal from the id.
			#gene_id = values[0].split(" ")
			#gene_id = gene_id[2]

			# Split the reads value, to separate the literal from the reads.
			#reads = values[1].split(" ")
			#reads = reads[2]
			reads = ""
			for i in counts:
				i = i.split(" ")
				i = i[2]
				reads += i + "\t"
			row += reads + "\n"
			
			#if self.gtf_index_count < 2:

				# Create row string.
				#row = k[1] + "\t" + gene_id + "\t" + count + "\t"

			# Write to the CSV file.
		self.transcriptExpressionCSV.write(row)

def main():
	p = ParseGTF("/groups/umcg-wijmenga/tmp04/umcg-eschutte/projects/gluten_specific_Tcells/output/FluxCapacitor/")
	p.GTFReader()
	p.assembleTranscriptExpressionMatrix()

if "__main__" == __name__:
	main()
