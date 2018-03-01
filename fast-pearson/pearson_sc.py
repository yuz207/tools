# Init
import numpy
import csv
import time
import os
import scipy
import itertools
from scipy.stats import pearsonr

# Calculate Pearson correlation pairwise over the rows (reg elements) of
# the imported signal arrays
start_time = time.time()

# Define output file location
# GIVE A LOCATION TO SAVE YOUR FILE ANTHONY
outfile_loc = 'pearson.txt'

# Define your input array of vectors
# MAKE CURR_SIGNAL = YOUR ARRAY HERE
curr_signal = numpy.random.random((10, 1000))
print curr_signal.shape
rows = curr_signal.shape[0]

# Init array to store pairwise correlation array
pearson_array = []

# Calculate Pearson row by row
means = curr_signal.mean(axis=1)[(slice(None, None, None), None)]
curr_signal_m = curr_signal - means
curr_signal_ss = numpy.sqrt(scipy.stats.ss(curr_signal_m, axis=1))

# Write rows to pearson.txt; also export each row independently as its own
# file, with row index # as filename
with open(outfile_loc, 'wb', 50000000) as pearson_out:
    csv_writer = csv.writer(pearson_out, delimiter="\t")
    for j in xrange(rows):
        # Calculate current row of Pearson correlations, p_row, to be written
        temp = numpy.dot(curr_signal_m, curr_signal_m[j].T)
        p_row = temp / (curr_signal_ss * curr_signal_ss[j])
        # print j, len(p_row)
        #p_row = numpy.around(p_row, decimals=10)
        # pearson_array.append(p_row)
        #p_row = [str(val) for val in p_row]
        # csv_writer.writerow(p_row)
        p_row = numpy.array(p_row)
        p_row = numpy.reshape(p_row, (1, -1))
        numpy.savetxt(pearson_out, p_row, fmt='%10.8f', delimiter='\t')
        #numpy.save(pearson_out, p_row)

print ("--- %s seconds ---" % (time.time() - start_time))
