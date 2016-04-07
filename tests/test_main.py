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
# SRR534843_100X_31mers.hgram
# SRR403105_both.fastq.31mers.hist
# SRR2880607_both_new100X.fastq.31mers.hist
# SRR1649436_100X_both.fastq.31mers.hist
# SRR1009262_both_100X.fastq.31mers.hist
# ERR560569_100X_both.fastq.31mers.hist
# yeast_sim.bfast.fastq.31mers.hist

class Test_find_start_repeat_kmers(unittest.TestCase):

	def run_find_start_repeat_kmers_test(self, hist_dict, lower_cutoff, upper_cutoff):
		estimated_start = grin_main.find_start_repeat_kmers(hist_dict)
		self.assertTrue(lower_cutoff <= estimated_start <= upper_cutoff,
				msg = "Failed. Start incorrectly estimated as %d" %(estimated_start))
	
	def test_fake1(self):
		with open("fake1_61mers.hgram", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 71, 81)

	def test_fake2(self):
		with open("fake2_61mers.hgram", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 71, 81)

	def test_celegans(self):
		with open("celegans_61mers.hgram", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 27, 37)

	def test_yeast(self):
		with open("yeast_51mers.hgram", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 64, 74)
	
	def test_SRR534843_100X_31mers(self):
		with open("SRR534843_100X_31mers.hgram", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 87, 97)

	def test_SRR403105_both(self):
		with open("SRR403105_both.fastq.31mers.hist", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 65, 75)

	def test_SRR2880607_both_new100X(self):
		with open("SRR2880607_both_new100X.fastq.31mers.hist", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 104, 114)

	def test_SRR1649436_100X_both(self):
		with open("SRR1649436_100X_both.fastq.31mers.hist", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 90, 100)

	def test_SRR1009262_both_100X(self):
		with open("SRR1009262_both_100X.fastq.31mers.hist", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 87, 97)

	def test_ERR560569_100X_both(self):
		with open("ERR560569_100X_both.fastq.31mers.hist", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 103, 113)

	def test_yeast_sim(self):
		with open("yeast_sim.bfast.fastq.31mers.hist", 'r') as f:
			hist_dict = grin_main.create_hist_dict(f)
			Test_find_start_repeat_kmers.run_find_start_repeat_kmers_test(self, hist_dict, 61, 71)
	
if __name__ == "__main__":
	unittest.main()
