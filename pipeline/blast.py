import subprocess as proc
import pandas as pd
import numpy as np

def create_blastdb(fasta_path, title, out_path,
                   dbtype='nucl', executable_path='makeblastdb'):
    """Creates a blast database by calling the makeblastdb executable.

    Parameters
    ----------
    fasta_path : str
    title : str
    out_path : str
    dbtype: str, optional

    Returns
    -------
    subprocess.CompletedProcess
    """
    blast_cmd = '{exec_path} -in {fasta_path} ' \
                '-input_type fasta -dbtype {dbtype} -title {title} ' \
                '-out {out_path}'
    cmd = blast_cmd.format(exec_path=executable_path,
                           fasta_path=fasta_path,
                           dbtype=dbtype,
                           title=title,
                           out_path=out_path)
    return proc.run(cmd, shell=True, stdout=proc.PIPE)

def call_blast(task, query_path, db_path, out_path,
               outfmt='6 std qlen slen qcovs sstrand',
               max_seq=1, threads=6, executable_path='blastn'):
    """Runs a blastn search.

    Parameters
    ----------
    task : str
    query_path : str
    db_path : str
    out_path : str
    outfmt : str, optional
    max_seq : int, optional
    threads : int, optional

    Returns
    -------
    subprocess.CompletedProcess
    """
    blast_cmd = '{exec_path} -task {task} -query {query} ' \
                '-db {db} -out {out} -outfmt {outfmt} ' \
                '-max_target_seqs {max_seq} -num_threads {threads}'
    cmd = blast_cmd.format(exec_path=executable_path,
                           task=task,
                           query=query_path,
                           db=db_path,
                           out=out_path,
                           outfmt='"'+outfmt+'"',
                           max_seq=max_seq,
                           threads=threads)
    return proc.run(cmd, shell=True, stdout=proc.PIPE)

def blast_to_df(path,
                filter_by_eval=True, eval_threshold=1e-10, sep='\t',
                col_labels=('qaccver', 'saccver', 'pident', 'length',
                            'mismatch', 'gapopen', 'qstart', 'qend',
                            'sstart', 'send', 'evalue', 'bitscore',
                            'qlen', 'slen', 'qcovs', 'sstrand')):
    """Converts a tab-delimited blastn result to a pandas DataFrame.

    Parameters
    ----------
    path : str
    filter_by_eval : bool, optional
    eval_threshold : float, optional
    sep : str, optional
    col_labels : iterable of str, optional

    Returns
    -------
    pandas.DataFrame
    """
    df = pd.DataFrame.from_csv(path, sep=sep, header=None, index_col=None)
    df.columns = col_labels
    if filter_by_eval:
        return df[df['evalue'] < eval_threshold]
    return df

def summarize(group):
    """Reduces a group of entries into a single row.

    Parameters
    ----------
    group : pandas.DataFrame

    Returns
    -------
    pandas.Series
    """
    return pd.Series({
        'qaccver': np.max(group['qaccver']),
        'saccver': np.max(group['saccver']),
        'pident': np.sum(group['pident'] * \
                  (group['length'] / group['length'].sum())),
        'length': group['length'].sum(),
        'mismatch': group['mismatch'].sum(),
        'gapopen': group['gapopen'].sum(),
        'bitscore': group['bitscore'].sum(),
        'qlen': np.max(group['qlen']),
        'slen': np.max(group['slen']),
    })

def group_blastdf(df,
                  summarize_data=True, groupby=('qaccver', 'saccver'),
                  summary_col_labels=('qaccver', 'saccver', 'pident',
                                      'length', 'mismatch', 'gapopen',
                                      'bitscore', 'qlen', 'slen')):
    """Groups blastn output based on the query and subject ID's.

    Parameters
    ----------
    df : pandas.DataFrame
    groupby : iterable of str
    summary_col_labels : iterable of str

    Returns
    -------
    pandas.DataFrame
    """
    grp_df = df.groupby(groupby)[summary_col_labels]
    if summarize_data:
        return grp_df.apply(summarize)[summary_col_labels[2:]]
    return grp_df
