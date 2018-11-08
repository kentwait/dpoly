OUTFMT = '"6 std qlen slen qcovs sstrand"'
MAX_SEQ = 1
THREAD_COUNT = 6
BLAST_CMD = '! blastn -task {task} -query {query} ' \
            '-db {db} -out {out} -outfmt {outfmt} ' \
            '-max_target_seqs {max_seq} -num_threads {threads}'
