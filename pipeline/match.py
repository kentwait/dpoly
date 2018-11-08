import subprocess as proc
from Bio import SeqIO
import pandas as pd

def match_by_fbtr(fbtr_list1, fbtr_list2):
    """Match FlyBase transcript ID's between 2 lists.

    Parameters
    ----------
    fbtr_list1 : list
    fbtr_list2 : list

    Returns
    -------
    list
        FlyBase transcript ID's found in both lists
    """
    fbtr_matches = set()
    for k in fbtr_list1:
        if k in fbtr_list2:
            fbtr_matches.add(k)
    return list(fbtr_matches)

def reciprocal_match_by_fbtr(fbtr_list1, fbtr_list2):
    """Reciprocally match FlyBase transcript ID's between 2 lists.

    Parameters
    ----------
    fbtr_list1 : list
    fbtr_list2 : list

    Returns
    -------
    list
        FlyBase transcript ID's found in both lists in both
        forward and reverse matching.
    """
    forward_match = set(match_by_fbtr(fbtr_list1, fbtr_list2))
    reverse_match = set(match_by_fbtr(fbtr_list2, fbtr_list1))
    return forward_match.intersection(reverse_match)

def call_blast(task, query_path, db_path, out_path,
               outfmt='6 std qlen slen qcovs sstrand',
               max_seq=1,
               threads=6):
    """

    Parameters
    ----------
    task : str
    query_path : str
    db_path : str
    out_path : str
    outfmt : str
    max_seq : int
    threads : int

    Returns
    -------
    int
    """
    blast_cmd = 'blastn -task {task} -query {query} ' \
                '-db {db} -out {out} -outfmt {outfmt} ' \
                '-max_target_seqs {max_seq} -num_threads {threads}'
    cmd = blast_cmd.format(task=task,
                           query=query_path,
                           db=db_path,
                           out=out_path,
                           outfmt='"'+outfmt+'"',
                           max_seq=max_seq,
                           threads=threads)
    return proc.run(cmd, shell=True, stdout=proc.PIPE)

def blast_to_df(path,
                filter_by_eval=True,
                eval_threshold=1e-10,
                sep='\t',
                col_labels=['qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'qlen', 'slen', 'qcovs', 'sstrand']):
    """Converts a tab-delimited blastn result to a pandas DataFrame.

    Parameters
    ----------
    path : str
    filter_by_eval : bool
    eval_threshold : float
    sep : str
    col_labels : list of str

    Returns
    -------
    pandas.DataFrame
    """
    df = pd.DataFrame.from_csv(path, sep=sep, header=None, index_col=None)
    df.columns = col_labels
    if filter_by_eval:
        return df[df['evalue'] < eval_threshold]
    return df