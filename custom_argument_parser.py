import sys
import argparse


class CustomParser(argparse.ArgumentParser):

	def generate_usage_message(self):
		return "usage: grin [-h] [-v] [-a] [-i] [-c cutoff [cutoff ...]] -f file [file ...]\n\n"

	def generate_help_message(self):
		store_str = self.generate_usage_message()
		store_str += "required arguments:\n"
		store_str += "-f, --file : one or more input histogram files\n\n"

		store_str += "optional arguments:\n"
		store_str += "-v, --verbose : print more output\n"
		store_str += "-c, --cutoffs : list of manual cutoffs for start of repetitive k-mers\n"
		store_str += "-a, --analyzer : run kmerspectrumanalyzer and compute GRI from that histogram instead\n"
		store_str += "-i, --ignore-errors : estimate the erroneous k-mers and do not include them in the total k-mer count\n"
		store_str += "-h, --help : print this message\n\n"

		store_str += "notes:\n"
		store_str += "* if specified, the number of manually specified cutoffs must equal the "
		store_str += "number of input files\n\n"

		return store_str

	def print_usage(self, file = None):
		if file is None:
			file = sys.stdout
		self._print_message(self.generate_usage_message(), file)

	def print_help(self, file = None):
		if file is None:
			file = sys.stdout
		self._print_message(self.generate_help_message(), file)
