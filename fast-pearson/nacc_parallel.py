# Init
import numpy, csv, time, os, itertools, scipy.stats

# Function to get the number of lines in a .bed file (aka number of feature coordinates)
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

## Define directories and paths of reference files and outputs
proj_dir_name = raw_input('Project directory name: ');
root = '../..';

# Define species assemblies to compare and classes of regulatory elements
species = ['hg19','mm9'];
#reg = 'RNA-seq2';
reg = raw_input('reg type: ');

# Define output directories, create them if they do not exist: out_dir/asm/reg_element
out_dir = root+'/nacc_analysis/'+proj_dir_name;


from itertools import izip, repeat, starmap
from multiprocessing import Pool, Queue, Process, Manager, cpu_count
print 'num_cpus: ' + str(cpu_count())

##
def par_pearson(j, curr_signal):
    # print j, curr_signal.shape

    # Current (partial) row of Pearson correlations, rs, to be written, containing values of upper triangular matrix
    #temp = numpy.dot(curr_signal_m[j:],curr_signal_m[j].T);
    #rs = temp / (curr_signal_ss[j:]*curr_signal_ss[j]);
    # Calculate current row of Pearson correlations, rs, to be written
    temp = numpy.dot(curr_signal_m[:], curr_signal_m[j].T)
    rs = temp / (curr_signal_ss[:] * curr_signal_ss[j])
    # print j, len(rs)
    #rs = numpy.around(rs, decimals=10)
    rs = [str(val) for val in rs]
    # csv_writer.writerow(rs)

    # If desired, write each line to own file in pearson subfolder
    with open(curr_array_folder + '/pearson/' + str(j) + '.txt', 'wb') as upper_row:
        row_writer = csv.writer(upper_row, delimiter="\t")
        row_writer.writerow(rs)

    return None

# Pool.starmap does not exist in Py 2.7; use itertools.starmap instead
def par_pearson_star(in_args):
    return par_pearson(*in_args)

##
for i in range(2):
    start_time = time.time()

    curr_asm = species[i]
    curr_array_folder = '/'.join([out_dir, reg, curr_asm, 'array'])
    if not os.path.exists(curr_array_folder + '/pearson'):
        os.makedirs(curr_array_folder + '/pearson')

    # Get signal value array from all_signal
    curr_signal = numpy.genfromtxt(
        '/'.join([out_dir, reg, curr_asm, 'array', 'avg_signal_norm.txt']))
    print curr_signal.shape
    #curr_signal = numpy.random.random((5000,5000))
    rows = curr_signal.shape[0]

    # Calculate Pearson row by row
    ms = curr_signal.mean(axis=1)[(slice(None, None, None), None)]
    curr_signal_m = curr_signal - ms
    # print curr_signal_m
    curr_signal_ss = numpy.sqrt(scipy.stats.ss(curr_signal_m, axis=1))

    # Define arguments for parallel_pearson
    job_args = zip(range(0, rows), itertools.repeat(curr_signal))

    if __name__ == '__main__':
        pool = Pool(processes=6)
        #pool.starmap(parallel_pearson, job_args)
        pool.map(par_pearson_star, job_args)
        # print len(all_partial_pearsons), all_partial_pearsons[0].shape
        pool.close()

    print ("--- %s seconds ---" % (time.time() - start_time))

