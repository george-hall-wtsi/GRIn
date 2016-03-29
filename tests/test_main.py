import unittest
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
import grin as grin_main

# Test Files:
# celegans_61mers.hgram fake1_61mers.hgram yeast_51mers.hgram fake2_61mers.hgram 
# SRR534843_100X_15mers.hgram SRR534843_100X_31mers.hgram SRR534843_100X_51mers.hgram 
# SRR534843_100X_71mers.hgram

class Test_find_start_repeat_kmers(unittest.TestCase):

	def run_find_start_repeat_kmers_test(self, hist_dict, lower_cutoff, upper_cutoff):
		estimated_start = grin_main.find_start_repeat_kmers(hist_dict)
		self.assertTrue(lower_cutoff <= estimated_start <= upper_cutoff,
				msg = "Failed. Start incorrectly estimated as %d" %(estimated_start))
	
	def test_fake1(self):
		with open("fake1_61mers.hgram", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 73, 79)

	def test_fake2(self):
		with open("fake2_61mers.hgram", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 73, 79)

	def test_celegans(self):
		with open("celegans_61mers.hgram", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 27, 33)
	
	def test_yeast(self):
		with open("yeast_51mers.hgram", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 66, 72)

	def test_SRR534843_100X_15mers(self):
		with open("SRR534843_100X_15mers.hgram", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 104, 110)

	def test_SRR534843_100X_31mers(self):
		with open("SRR534843_100X_31mers.hgram", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 90, 96)

	def test_SRR534843_100X_51mers(self):
		with open("SRR534843_100X_51mers.hgram", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 67, 73)

	def test_SRR534843_100X_71mers(self):
		with open("SRR534843_100X_71mers.hgram", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 40, 46)

if __name__ == "__main__":
	unittest.main()
