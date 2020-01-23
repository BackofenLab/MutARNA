#!/usr/bin/env python2
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import os
import pandas as pd
from tempfile import NamedTemporaryFile
from subprocess import Popen, PIPE
from io import StringIO
import time
from math import ceil, floor
import re


# rase_root_dir = os.path.join(os.environ['HOME'], 'repositories/RaSE/')
# rase_src_dir = os.path.join(rase_root_dir, 'code')
# sys.path = [rase_src_dir] + sys.path

# eden_root_dir = os.path.join(os.environ['HOME'], 'repositories/EDeN/')
# eden_src_dir = os.path.join(eden_root_dir)
# sys.path = [eden_src_dir] + sys.path


# edenrna_root_dir = os.path.join(os.environ['HOME'], 'repositories/eden_rna/')
# edenrna_src_dir = os.path.join(edenrna_root_dir)
# sys.path = [edenrna_src_dir] + sys.path


def main_SNV_wrapper(argv):
    print (sys.argv)
    input_file = sys.argv[1]
    output_file_prefix = sys.argv[2]
    window = int(sys.argv[3])
    if len(sys.argv) > 4:
        num_splits = int(sys.argv[4])
        split_id = int(sys.argv[5])
    else:
        num_splits = 1
        split_id = 0

    start_time = time.time()

    #initiate empty dataframes for RNAsnp
    df_rnasnp = pd.DataFrame(columns=['SNP', 'w', 'Slen', 'GC', 'interval', 'd_max', 'p-value', 'interval.1', 'r_min', 'p-value.1', 'ID'])

    #initiate empty dataframes for remuRNA
    df_remurna = pd.DataFrame(columns=['SNP', 'MFE(wt)', 'MFE(mu)', 'dMFE', 'H(wt||mu)', 'GCratio', 'ID'])

    rase_scores = []
    lcount = 0
    fasta_sequences = list(SeqIO.parse(open(input_file), 'fasta'))
    total_size = len(fasta_sequences)
    ranges =  list(range(0, total_size, int(ceil(total_size/float(num_splits)))))
    print('runner.py args:' + ' '.join(sys.argv))
    print('rangesA: ', ranges)
    if ranges[-1] != total_size:
        ranges.append(total_size)
    print ('rangesB: ', ranges)
    print ('runner on range: ', ranges[split_id], ranges[split_id+1])
    runRnasnp = True#False
    runRemurna = True#False
    runRase = False
    for fasta in fasta_sequences[ranges[split_id]: ranges[split_id+1]]:
        lcount += 1
        print ('\r{}..' .format(lcount),) 
        id, desc, sequence = fasta.id,fasta.description, str(fasta.seq)
        #extract snp from description or id
        match = re.search("(\w\d+\w)$", desc)
        if match is None:
            raise RuntimeError('SNP tag not found for desc:{}'.format(desc))

        snp = [match.group(1)]
        print (snp,)
        tmp_seq_fa = NamedTemporaryFile(suffix='.fa', delete=False)
        tmp_seq_fa.write(">" +desc + "\n")
        tmp_seq_fa.write(sequence)
        tmp_seq_fa.close()

        #run RNAsnp
        if runRnasnp:
            res_rnasnp = run_RNAsnp(tmp_seq_fa.name,snp,window)
            res_rnasnp = res_rnasnp.assign(ID=id)
            df_rnasnp = df_rnasnp.append(res_rnasnp, ignore_index=True)

        #run remuRNA
        if runRemurna:
            res_remurna = run_remuRNA(tmp_seq_fa.name,snp, window)
            res_remurna = res_remurna.assign(ID=id)
            df_remurna = df_remurna.append(res_remurna, ignore_index=True)

        ##run RaSE
        if runRase:
            rase_scores = rase_scores+ [[id,snp[0],run_RaSE(sequence,snp, window=window)]]
                    
        # remove temp file
        os.remove(tmp_seq_fa.name)

    output_file_prefix += "_"+str(split_id)
    if runRnasnp:
        df1 = df_rnasnp.set_index('ID')
        df1['tool-parameters:window={}'.format(window)] = ''
        df1.to_csv(path_or_buf=output_file_prefix+"_rnasnp.csv", sep="\t")

    if runRemurna:
        df2 = df_remurna.set_index('ID')
        df2['tool-parameters:window={}'.format(window)] = ''
        df2.to_csv(path_or_buf=output_file_prefix+"_remurna.csv", sep="\t")

    if runRase:
        df_rase= pd.DataFrame(rase_scores, columns=['ID','SNP','Score'])
        df_rase =df_rase.set_index('ID')
        df_rase['tool-parameters:window={}|avg_bp_prob_cutoff=0.01|hard_threshold=0.5|max_num_edges=3'.format(window)]='' 
        df_rase.to_csv(path_or_buf=output_file_prefix+"_rase.csv",sep="\t")
    

    print("--- %s seconds ---" % (time.time() - start_time))

def run_remuRNA(wild_fa, snp_tags, window=None):
    """
    A python wrapper invoking remuRNA tool.
    Please check remuRNA documentation for details.
    Call example: run_remuRNA('./wild.fa', ['G20C'])
    Parameters
    ----------
    wild_fa : str
        Fasta file name containing one RNA sequence
    snp_tags : list
        Set of SNP tags required to be evaluatued on the input sequence.
        Warning: remuRNA accepts only a single tag in each call.

    Returns
    -------
    dataframe
        Pandas table of standard output

    """
    assert(len(snp_tags) == 1)
    if not os.path.isfile(wild_fa):
        raise RuntimeError("Input fasta %s does not exist" % in_fasta_file)

    # Write RNA sequence and SNP tags to a temporary file, 
    tmp_fa = NamedTemporaryFile(suffix='.fa', delete=False)
    with open(tmp_fa.name, 'w') as out_fa_handle:
        with open(wild_fa) as in_fa_handle:
            for line in in_fa_handle.readlines():
                out_fa_handle.write(line)
        out_fa_handle.write('\n'.join(['*'+tag for tag in snp_tags]))

    
    # Make a shell command line
    cmd = '$(which remuRNA) {}  '.format(tmp_fa.name)
    params = '-p=1 '
    if window is not None:
        params += '-w={}'.format(int(window))
    cmd += params
    # print cmd
    p = Popen(cmd , stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    if err:
        print("Error in calling remuRNA\n{}\n{}\n".format(out, err))
        df = pd.DataFrame({'SNP':{0:snp_tags[0]}, 'error':{0:'True:{}'.format(err.decode('utf-8'))}, })
    else:
        df = pd.read_table(StringIO(out.decode("utf-8")))
    df['remurna_params'] = params
    return  df
    #return out

def run_RNAsnp(wild_fa, snp_tags, window=None, plfold_W=None, plfold_L=None,  mode=None, param_string=''):
    """
    A python wrapper invoking RNAsnp tool.
    Please check RNAsnp documentation for details.

    Call example: run_RNAsnp('./wild.fa', ['G20C'])
    Parameters
    ----------
    wild_fa : str
        Fasta file name containing one RNA sequence
    snp_tags : list
        Set of SNP tags required to be evaluatued on the input sequence
    window : int
        Length of flanking sequence on either side of. If None, the RNAsnp window (-W) size. Windows larger than 800 will be passed as 800.

    Returns
    -------
    dataframe
        Pandas table of standard output

    """

    # Write SNP tags to a temporary file, 

    snp_file = NamedTemporaryFile(delete=False)

    with open(snp_file.name, 'w') as out_handle:
        out_handle.write('\n'.join(snp_tags))
        
    if not os.path.isfile(wild_fa):
        raise RuntimeError ("Input fasta %s does not exist" % wild_fa)
    

    # Make a shell command line
    cmd = 'RNAsnp -f {} -s {} '.format(wild_fa, snp_file.name)
    params = ""
    if mode is not None:
        assert mode in [1, 2, 3]
        params += '-m {} '.format(mode) 
    if plfold_W is not None:
        params  += '-W {} '.format(plfold_W)
    if plfold_L is not None:
        params  += '-L {} '.format(plfold_L)

    if window is not None:
        if window > 800:
            print ("WARNING RNAsnp window reduced to max possible: 800")
            window = 800
        params += '--winsizeFold {} '.format(int(window))
    cmd += params + param_string
    #print cmd 
    p = Popen(cmd, stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    if err:
        #raise RuntimeError("Error in calling RNAsnp\n{}\n{}\n".format(out, err))
        print ("Error in calling RNAsnp\n{}\n{}\n".format(out, err))
    #print out

    # os.remove(snp_file.name)
    print(params)
    out_cleaned = ""
    for line in out.decode('utf-8').split('\n'):
        line = line.strip()
        if len(line)==0:
            continue
        if 'error' in line.lower():
            print("RNASNP returned error for: {} message is:{}".format(wild_fa, line))
        elif 'warning' in line.lower():
            print("ERROR: RNASNP returned warning for: {} message is:{}".format(wild_fa, line))
        else:
            out_cleaned += line+"\n"
    if len(out_cleaned) != 0:
        print(out_cleaned)
        df_RNAsnp = pd.read_table(StringIO(out_cleaned))
    else:
        df_RNAsnp = pd.DataFrame({'SNP':{0:snp_tags[0]}, 'error':{0:'True:{}'.format(err.decode('utf-8'))}, })
    df_RNAsnp['rnasnp_params'] = params + param_string
    return  df_RNAsnp #.add_suffix('RNAsnp:')
    #return out_cleaned



if __name__ == "__main__":
   main_SNV_wrapper(sys.argv[1:])
