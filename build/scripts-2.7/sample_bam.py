#!/usr/local/anaconda2/bin/python

# coding: utf-8

# In[ ]:

#__version__= 0.3.0


# In[1]:

import argparse, pysam, logging, subprocess,  tempfile, re, os


# In[6]:

def parse_args():
    parser = argparse.ArgumentParser(description='''
                                     This script takes a bam filename and a series of requested read numbers.
                                     It then samples mapped read pairs and generates new subsampled bam files
                                     with the desired number of reads. If mapped reads pairs are present (and
                                     uniquely mapped), the pair is randomly selected and written.''')
    parser.add_argument('bam_filename', help='Name of bam file to subsample')
    parser.add_argument('-o', '--output_dir', help='Optional output directory. If omitted, bams are written to the current directory')
    parser.add_argument('sample_x_n', nargs='+', help='''
    Number of reads to sample and number of sampled bams to generate. Multiple entries can be given.
    For example, sample_bams.py test.bam 1000x5 10000x2 will generate 5 bams with 1000 reads and 2 bams with 10000 reads
    ''')
    args = parser.parse_args()

    args.output_dir = '' if args.output_dir is None else args.output_dir
    print "FOO!"
    return args
    #return dict(args)


# In[3]:

def sample_and_write_file(bam_in_filename, bam_out_filename, sample_probability):
    """
    This function randomly selects read pairs based on the provided sample probability and writes these to bam_out_filename.
    Every time it encounters a read thast is not in encountered_readnames_set, it selects it based on the provided probability
    and then adds its query name to the encountered_readnames_set if selected. Whenever a read is encountered it is automatically
    written if it's query name is in encountered_readnames_set, which is followed by removal of the read query name from
    encountered_readnames_set.

    Potential problems to be aware of:
    1) Query names are matched exactly, so if there are read direction suffixes (e.g. /1 and /2), then no mates will ever be found.
       This will result in an ever expanding encountered_readnames_set --possibly resulting in exhausting the system's memory-- and
       and a set of reads with some paired (by chance selection) and some not paired.
    2) If read pairs are mapped more non-uniquely, then each pair of encountered read names will be treated as unique paired reads.
       This is because after encountering the matched query, the query name will be cleared from encountered_readnames_set so the
       next encounter of a read from this pair will be treated as the first read in a newly enountered read pair.

    """
    import random
    encountered_readnames_set = set()
    bam_in = pysam.AlignmentFile(bam_in_filename, 'rb')
    bam_out = pysam.AlignmentFile(bam_out_filename, 'wb', template = bam_in)

    for _read in bam_in.fetch():
        if _read in encountered_readnames_set:
            bam_out.write(_read)
            #Remove from set ONLY if read pair occurs ONCE in bam, so not if  using gmapper for example
            encountered_readnames_set.remove(_read.qname)
        else:
            if random.uniform(0,1) <= sample_probability:
                bam_out.write(_read)
                encountered_readnames_set.add(_read.qname)
    return True


# In[17]:

def get_sorted_list_of_sample_bam_lengths_list(sample_x_n_list):
    '''
    This function takes a list of sample_x_n (e.g. [1000x2, 100x3]) and converts
    it to a sorted list of lists of desired bam read lengths.
    For the above example, the returned list would be [[100, 100, 100], [10000, 10000]]
    This function has no side effects
    '''
    lengths_list = []
    for _element in sample_x_n_list:
        nreads_to_nsamples_list = map(int, _element.split("x"))
        lengths_list.append([nreads_to_nsamples_list[0]]*nreads_to_nsamples_list[1])
    return sorted(lengths_list)


# In[21]:

def get_mapped_paired_read_count(bam_filename):
    """
    This function determines the number of mapped reads in a bam file.
    It calls samtools flagstat, reads stdout, then pulls the first word (the read count) from the third line after verifying it contains the
    word 'mapped'. The count is returned as an integer
    """
    flagstat_list = subprocess.check_output(['samtools', 'flagstat', bam_filename]).split("\n")
    if "supplementary" in flagstat_list[2]:
        assert re.search('mapped', flagstat_list[4])
        return int(flagstat_list[4].split(" ")[0])
    else:
        assert re.search('mapped', flagstat_list[2])
        return int(flagstat_list[2].split(" ")[0])


# In[85]:


def sort_and_index_bamfile(bam_filename):
    subprocess.check_call(['samtools', 'sort', '{}'.format(bam_filename), '{}'.format(re.sub('\..*$', '', ))])
    subprocess.check_call(['samtools', 'index', '{}'.format(bam_filename)])

if __name__ == '__main__':
    args = parse_args()
    print args
    print 'hello!'
    n_mapped_reads = get_mapped_paired_read_count(args.bam_filename)
    sorted_bam_lengths_list_list = get_sorted_list_of_sample_bam_lengths_list(args.sample_x_n)
    for same_size_list in sorted_bam_lengths_list_list:
        num_samples_of_same_length = len(same_size_list)
        assert(len(set(same_size_list)) == 1)
        bam_length = same_size_list[0]
        if bam_length > n_mapped_reads:
            print('The number of reads to sample must be less than the number of mapped reads in the bam [n={}]. Not generating file'.format(str(n_mapped_reads)))
            continue
        for _bam_length_index in range(0, num_samples_of_same_length):
            bam_out_filename = os.path.join(args.output_dir,re.sub(".bam", "", os.path.basename(args.bam_filename)) + "_sample_" + "{:0>5d}".format(_bam_length_index) + "_n" + str(bam_length) + ".bam")
            sample_and_write_file(args.bam_filename, bam_out_filename, float(bam_length)/(n_mapped_reads))
            sort_and_index_bamfile(bam_out_filename)

