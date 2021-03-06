#!/bin/env python

import argparse, pysam, logging, subprocess,  tempfile, re, os, sys


# In[6]:

def parse_args():
    parser = argparse.ArgumentParser(description='''
                                     This script subsamples a series of fastq
                                     files and outputs the subsampled versions.
                                     The sampling is done such that the indices
                                     in the original file are preserved. For
                                     example, entry 3 in the subsampled files
                                     will correspond the nth entries from the
                                     original files.
                                     ''')
    parser.add_argument('-o', '--output_dir', help='''Optional output directory.
                        If omitted, sampled files are written to the current
                        directory''')
    parser.add_argument('n_output_reads', help='''The approximate number of reads to
                        output''')
    parser.add_argument('n_input_reads', help='''The approximate number
                        of reads in the input fastq file(s)''')
    parser.add_argument('fastq_filenames', nargs='+', help='''
                        all fastq filenames to sample from. These will be
                        grouped such that the same indices will be selected
                        from each file during the sampling process.
    ''')
    args = parser.parse_args()

    args.output_dir = '' if args.output_dir is None else args.output_dir
    return args

def read_sequence(file_obj):
    '''
    This function takes a file object as input, reads 4 lines of input, and
    returns the concatenated lines. The 3rd line is checked to make sure it
    begains with a '+' character, as is usual with the Fastq format. If this
    ceck fails, the program crashes
    '''

    seq_lines = [file_obj.readline() ,file_obj.readline() ,file_obj.readline(),
                file_obj.readline()]
    print(seq_lines)
    if not re.search('^\+', seq_lines[2]):
        sys.exit('''
                 Are you sure this is a fastq file? Where's the '+' character?''')
    else:
        return ''.join(seq_lines)

def sample_and_write_files(fastq_filenames_list, output_directory, sample_probability):
    """
    This function iterates over each of the fastq_files in
    fastq_filenames_list, and samples reads according to sample_probability.
    Sampled reads are output to the output_directory. The sampled file names
    are modified such that .fastq is converted to _sub.fastq
    """
    import random
    file_obj_list = []
    for _filename in fastq_filenames_list:
        output_filename = os.path.join(output_directory, re.sub(".fastq",
                                                                "_sub.fastq",
                                                                os.path.basename(_filename)))
        file_obj_list.append([open(_filename, 'r'), open(output_filename, 'w')])
    
    FILE_END  = False
    while not FILE_END:
        select_bool = True if (random.uniform(0,1) <= sample_probability) else False
        for _file_objs in file_obj_list:
            try:
                seq_line = read_sequence(_file_objs[0])
                if select_bool:
                    _file_objs[1].write(seq_line)
            except IOError as e:
                FILE_END = True



if __name__ == '__main__':
    args = parse_args()
    print(args)
    sample_probability = float(args.n_output_reads)/float(args.n_input_reads)
    sample_and_write_files(args.fastq_filenames, args.output_dir,
                           sample_probability)
    print("DONE!")
