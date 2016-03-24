import sys
import scipy.signal
import numpy as np

def find_start_repeat_kmers(hist_dict):

	window_size = 5
	window = np.ones(int(window_size))/float(window_size)
	moving_average = np.convolve(hist_dict.values(), window, 'same')
	smoothed_data = dict(zip(hist_dict.keys(), [int(x) for x in moving_average]))

	order_num = 2

	min_list = scipy.signal.argrelextrema(np.array(smoothed_data.values()), np.less_equal, 
		order = order_num)[0].tolist()
	max_list = scipy.signal.argrelextrema(np.array(smoothed_data.values()), np.greater_equal, 
		order = order_num)[0].tolist()

	#print "min_list:" , min_list
	#print "max_list:" , max_list

	return (2 * max_list[1]) - min_list[1]

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
		start = find_start_repeat_kmers(hist_dict)[2]

		print start

		total_number_kmers = sum((a * b) for (a, b) in hist_dict.items())
		print total_number_kmers

		number_repetitive_kmers = 0
		for (a, b) in hist_dict.items():
			if (a >= start):
				number_repetitive_kmers += (a * b)

		print number_repetitive_kmers
