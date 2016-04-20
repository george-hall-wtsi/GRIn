import unittest
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import grin as grin_main


# Test Files:
# celegans_61mers.hgram
# fake1_61mers.hgram
# yeast_51mers.hgram
# fake2_61mers.hgram


class Test_calculate_GRI(unittest.TestCase):
	
	# grin_main.calculate_gri(hist_dict, verbose, error_cutoff, upper_bound, 
	# 			start_repetitive_kmers = 0):                                 

	def __init__(self, *args, **kwargs):
		super(Test_calculate_GRI, self).__init__(*args, **kwargs)
		self.in_file = "celegans_61mers.hgram"

	def create_hist_dict(self, in_file):
		with open(in_file, 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
		
		return hist_dict

	def test_default(self):

		hist_dict = self.create_hist_dict(self.in_file)
		gri = grin_main.calculate_gri(hist_dict, False, 0, None)
		self.assertEqual(gri, 0.11398022014165018)

	def test_error_cutoff(self):

		hist_dict = self.create_hist_dict(self.in_file)
		gri = grin_main.calculate_gri(hist_dict, False, 4, None)
		self.assertEqual(gri, 0.12392366647125003)

	def test_upper_cutoff(self):

		hist_dict = self.create_hist_dict(self.in_file)
		gri = grin_main.calculate_gri(hist_dict, False, 0, 1000)
		self.assertEqual(gri, 0.08961701517325185)

	def test_both_cutoffs(self):

		hist_dict = self.create_hist_dict(self.in_file)
		gri = grin_main.calculate_gri(hist_dict, False, 4, 1000)
		self.assertEqual(gri, 0.09766934730858438) 


if __name__ == "__main__":
	unittest.main()
