import sys

if len(sys.argv) != 4:
	print "usage: %s <start of first peak> <mid point first peak> <file>" %(sys.argv[0])
	sys.exit()

# 'a' is section before mid point, 'b' is section after
kmers_in_a = 0
kmers_in_b = 0
number_a = 0
number_b = 0

start_peak_1 = int(sys.argv[1])
midpoint = int(sys.argv[2])

with open(sys.argv[3], 'r') as f:
	for line in f.readlines():
		splat = [int(x) for x in line.strip().split()]
		if start_peak_1 <= splat[0] <= midpoint:
			kmers_in_a += (splat[0] * splat[1])
                        number_a += splat[1]
		elif splat[0] > midpoint:
			kmers_in_b += (splat[0] * splat[1])
                        number_b += splat[1]

	print "k-mers before midpoint =" , kmers_in_a
	print "k-mers after midpoint =" , kmers_in_b
	print "GRI =" , ((1.0 * kmers_in_b) / kmers_in_a)
        print "Other =" , ((1.0 * number_b) / number_a)
