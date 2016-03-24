import unittest
import sys

import calc_gri as calc_gri_main

# Test Files: celegans_61mers.hgram fake1_61mers.hgram yeast_51mers.hgram fake2_61mers.hgram

class Test_find_start_repeat_kmers(unittest.TestCase):
	
	def test_fake1(self):

		with open("fake1_61mers.hgram", 'r') as f:
			hist_dict = calc_gri_main.create_hist_dict(f)
			self.assertTrue(73 <= calc_gri_main.find_start_repeat_kmers(hist_dict) <= 79)

	def test_fake2(self):

		with open("fake2_61mers.hgram", 'r') as f:
			hist_dict = calc_gri_main.create_hist_dict(f)
			self.assertTrue(70 <= calc_gri_main.find_start_repeat_kmers(hist_dict) <= 76)

	def test_celegans(self):

		with open("celegans_61mers.hgram", 'r') as f:
			hist_dict = calc_gri_main.create_hist_dict(f)
			self.assertTrue(23 <= calc_gri_main.find_start_repeat_kmers(hist_dict) <= 29)

	def test_yeast(self):

		with open("yeast_51mers.hgram", 'r') as f:
			hist_dict = calc_gri_main.create_hist_dict(f)
			self.assertTrue(64 <= calc_gri_main.find_start_repeat_kmers(hist_dict) <= 70)


if __name__ == "__main__":
	unittest.main()
