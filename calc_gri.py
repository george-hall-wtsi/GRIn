import sys
import scipy.signal
import numpy as np

def find_start_repeat_kmers(hist_dict):

	order_num = 1

	min_list = scipy.signal.argrelextrema(np.array(hist_dict.values()), np.less_equal, 
		order = order_num)[0].tolist()
	max_list = scipy.signal.argrelextrema(np.array(hist_dict.values()), np.greater_equal, 
		order = order_num)[0].tolist()

	min_list_minimum = min(min_list)

	for maximum in sorted(max_list):
		if (maximum > min_list_minimum) and (maximum > 10):
			first_peak = maximum
			break

	# Return the point 'x' such that the first peak is equidistant 
	# between the first minimum and x
	return ((2 * first_peak) - min_list_minimum)

def create_hist_dict(in_file):
	hist_dict = {}
	for line in in_file.readlines():
		splat = [int(x) for x in line.strip().split()]
		hist_dict[splat[0]] = splat[1]

	return hist_dict


if __name__ == "__main__":

	if len(sys.argv) != 2:
		print "usage: %s <file>" %(sys.argv[0])
		sys.exit()

	with open(sys.argv[1], 'r') as f:
		hist_dict = create_hist_dict(f)
		start_repetitive_kmers = find_start_repeat_kmers(hist_dict)

		print "Start of repetitive k-mers" , start_repetitive_kmers

		total_number_kmers = sum((a * b) for (a, b) in hist_dict.items())
		print "Total number of k-mers" , total_number_kmers

		number_repetitive_kmers = 0
		for (a, b) in hist_dict.items():
			if (a >= start_repetitive_kmers):
				number_repetitive_kmers += (a * b)

		print "Number of repetitive k-mers" , number_repetitive_kmers

