import sys
import argparse


class CustomParser(argparse.ArgumentParser):

	def generate_usage_message(self):
		store_str = "usage: grin [-h] [-v] [-a] [-u upper-bound] [-i] [-e error-cutoff "
		store_str += "[error-cutoff ...]] [-c cutoff [cutoff ...]] -f file [file ...]\n"

		return store_str

	def generate_help_message(self):
		store_str = "\n"
		store_str += self.generate_usage_message()
		store_str += "\nrequired arguments:\n"
		store_str += "\t-f, --file              one or more input histogram files\n\n"

		store_str += "optional arguments:\n"
		store_str += "\t-v, --verbose           print more output\n"
		store_str += "\t-a, --analyzer          run kmerspectrumanalyzer and compute GRI "
		store_str += "from that histogram instead\n" 
		store_str += "\t-u, --upper-bound       upper bound after which k-mers do "
		store_str += "not contribute to total or repetitive k-mer counts\n"
		store_str += "\t-i, --ignore-error      estimate the erroneous k-mers and do "
		store_str += "not include them in the total k-mer count\n"
		store_str += "\t-e, --error-cutoffs     manual list of the end of the error peak "
		store_str += "(one per file)\n"
		store_str += "\t-c, --cutoffs           list of manual cutoffs for start of "
		store_str += "repetitive k-mers\n\n"
		
		store_str += "other arguments:\n"
		store_str += "\t-h, --help              print this message\n\n"

		store_str += "notes:\n"
		store_str += "\t* if specified, the number of manually specified cutoffs must "
		store_str += "equal the number of input files\n"
		store_str += "\t* if specified, the number of manually specified error cutoffs must "
		store_str += "equal the number of input files\n"
		store_str += "\t* --error-cutoffs and --ignore-error are mutually exclusive\n"
		store_str += "\t* --upper-bound must be a positive integer\n\n"

		return store_str

	def print_usage(self, file = None):
		if file is None:
			file = sys.stdout
		self._print_message(self.generate_usage_message(), file)

	def print_help(self, file = None):
		if file is None:
			file = sys.stdout
		self._print_message(self.generate_help_message(), file)
