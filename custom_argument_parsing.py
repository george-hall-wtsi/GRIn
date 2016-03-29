import sys
import argparse


def create_parser():
	parser = CustomParser(description = "Calculate the Genome Repeat Index")
	parser.add_argument("-v", "--verbose", action = "store_true", 
			help = "print more output")
	parser.add_argument("-c", "--cutoffs", nargs = '+', 
			help = "set manual cutoff (if specified, number of manual cutoffs must be the same\
			as the number of files)")
	parser.add_argument("-f", "--file", type = str, nargs = '+', required = True, 
			help = "input file(s)")

	return parser


class CustomParser(argparse.ArgumentParser):

	def generate_usage_message(self):
		return "usage: grin [-h] [-v] [-c cutoff [cutoff ...]] -f file [file ...]\n\n"

	def generate_help_message(self):
		store_str = self.generate_usage_message()
		store_str += "required arguments:\n"
		store_str += "-f, --file : one or more input histogram files\n\n"

		store_str += "optional arguments:\n"
		store_str += "-v, --verbose : print more output\n"
		store_str += "-c, --cutoffs : list of manual cutoffs for start of repetitive k-mers\n"
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


def parser_main():
	parser = create_parser()
	args = parser.parse_args()

	if args.cutoffs:
		if len(args.file) != len(args.cutoffs):
			print "ERROR: Need to have the same number of manual cutoffs as files"
			sys.exit(1)

	return args

