import os
import pandas as pd
from subprocess import Popen, PIPE
import re
import argparse

RCHIE_BIN = '/home/milad/1workspace/forked_repositories/r-chie/rchie.R'


def is_valid_file(file_name):
    if os.path.isfile(file_name):
        return os.path.abspath(file_name)
    else:
        raise FileNotFoundError(os.path.abspath(file_name))
def is_valid_directory(dir_name):
    if os.path.isdir(dir_name):
        return os.path.abspath(dir_name)
    else:
        raise FileNotFoundError(os.path.abspath(dir_name))



def run_intaRNA_csv(query_fa, target_fa, query_id, target_id, out_csv_file=None, shape_file=None, out_suffix=''):

    # Make a shell command line
    #subopt_csv = 'intarna-subopts_{}-{}{}.csv'.format(target_id, query_id, out_suffix)
    INTARNA_BIN = '' #os.path.join('') 
    INTARNA_ARGS =  '--outMode C --outCsvCols "id1,start1,end1,id2,start2,end2,E,E_hybrid" '\
    '--outOverlap Q --seedMaxUP 2 --temperature 37 '\
    '-n 100 '
    # '--qAccW 200 --qAccL 150 --tAccW 200 --tAccL 150 '\
    

    cmd = INTARNA_BIN + 'IntaRNA --target={} --query={} '.format(target_fa, query_fa)
    if out_csv_file:
        cmd += '--out={} '.format(out_csv_file)

    params = '{} '.format(INTARNA_ARGS)
    
        
    if shape_file is not None:
        cmd += '--qShape {} '.format(shape_file)
        params += '--qShapeMethod "Z" '

    print (cmd, params)
    if out_csv_file is not None:
        with open(out_csv_file, 'w') as csvhandle:
            p = Popen(cmd + params , stdin=PIPE, shell=True, stdout=csvhandle, stderr=PIPE)
            out, err = p.communicate()
            if err: 
                raise RuntimeError("Error in calling IntARNA\n{}\n".format(err))

    else:
        p = Popen(cmd + params , stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)
        out, err = p.communicate()
        if err or b"ERROR" in out: 
            raise RuntimeError("Error in calling IntARNA\n{}\n{}\n".format(out, err))
    
    #print (out)
    #df_subopts = pd.read_csv(subopt_csv,sep=';')

    #df_subopts['intarna_params:{}'.format(params)] = ''
    return  out_csv_file#df_subopts#, df_heatmap, spotprob_csv

# query_file, target_file = '../tmp/input_query.fa', '../tmp/input_target.fa'
# run_intaRNA_csv(query_file, query_file, "Q", "T", out_csv_file='../tmp/intarnaout.csv')

def csv_to_hlx(intarna_csv, len_target, shift_len_query, hlx_scale=1, hybrid_energy=False):

    df = pd.read_csv(intarna_csv, delimiter=";")
    
    intarna_rename_dict = {"target_start_pos":"start1","sRNA_start_pos":"start2", 
                       "target_end_pos":"end1","sRNA_end_pos":"end2",
                       "hybrid_energy":"E_hybrid", "energy":"E"}

    # rename webserver output to cmd column heading format
    df.rename(intarna_rename_dict, inplace=True, axis=1)
#     print(df.co)
    if hybrid_energy:
        value_col = 'E_hybrid' #'hybrid_energy'
    else:
        value_col = 'E' #'energy'
    df["target_hyb_len"] = df["end1"]-df["start1"]+1
    df["query_hyb_len"] = df["end2"]-df["start2"]+1
#     df_hlx = df[["target_start_pos","sRNA_start_pos","sRNA_hyb_len",value_col]].rename({'target_start_pos':'I',"sRNA_start_pos":'j',"sRNA_hyb_len":"length",value_col:"value"},axis=1).sort_values('value')
    df_hlx = df[["start1","start2","query_hyb_len",value_col]].rename({'start1':'I',"start2":'j',"query_hyb_len":"length",value_col:"value"},axis=1).sort_values('value')    
    df_hlx['j'] += shift_len_query 
    hlx_file = intarna_csv 
    if hybrid_energy:
        hlx_file+= '.hybrid.hlx'
    else:
        hlx_file += '.hlx'
    str_hlx = df_hlx.to_csv(sep='\t',index=None)
    with open(hlx_file,'w') as hlxout:
        hlxout.write('#'+str(int(len_target*hlx_scale))+'\n'+str_hlx)
    print("Converted {} to {}".format(intarna_csv, hlx_file))
    return hlx_file

# !head '../results/intarna/2019/WT-C100-C1000-5281772/predictions.sorted.csv'
# WT_hlx = csv_to_hlx('../results/intarna/2019/WT-C100-C1000-5281772/predictions.sorted.csv',hybrid_energy=False, context_len=1000)
# ! head $WT_hlx


def plot_rchie(WT_hlx, MUT_hlx,  out_name, energy_large_low_mid_high, convert_png=False, convert_svg=False):

    
    # RCHIE_BIN =  '/home/milad/1workspace/MutaRNA/repos/r-chie/rchie.R'

    color_ranges = ['\-12.5,\-10,\-7.5,\-5,\-4',
                        '\-25,\-20,\-15,\-10,\-8',
                        '\-40,\-35,\-30,\-20,\-10',
                        '\-50,\-40,\-30,\-20,\-15,\-10,\-8'
                    
                       ]
    
    color_range = color_ranges[energy_large_low_mid_high]
    
    rchie_cmd = RCHIE_BIN + ' '
    rchie_cmd += WT_hlx  + r' --format1 helix  --palette1 YlOrBr --group1 {}'.format(len(color_range.split(','))) + ' --value1 "'+color_range+ '" --rule1 4 ' 
    rchie_cmd += MUT_hlx + r' --format2 helix  --palette2 YlOrBr --group2 {}'.format(len(color_range.split(','))) + ' --value2 "'+color_range+ '" --rule2 4 '
    rchie_cmd += r' --basecol --legend2 --decreasing1=0 --decreasing2=0 --minimum1=-100 --maximum1=-4 --minimum2=-100 --maximum2=-4 '
    
    rchie_cmd += ' --pdf --output {}.pdf'.format(out_name)
#     rchie_cmd_png += ' --output {}.png'.format(out_name)
    
    if convert_png:
        rchie_cmd += ' ; magick {}.pdf {} {}.png'.format(out_name, "-resize 1500x1000", out_name)
    if convert_svg:
        rchie_cmd += ' ; pdf2svg {}.pdf {}.svg '.format(out_name, out_name)
        no_logo = True
        if no_logo:
            rchie_cmd += " &&  sed \'s/www.e.rna.org/ /g\' -i.bak  {}.svg".format(out_name)

    from subprocess import Popen, PIPE
    p = Popen(rchie_cmd, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    if (err):
        print(rchie_cmd)
        print('Warning/error in Rchie:\n{}'.format(err.decode('utf-8'))) 
    print (out.decode('utf-8'))
    

def plot_all(WT_csv, MUT_csv,  out_dir, len_target, shift_len_query, convert_png=False, convert_svg=True):
    out_name = os.path.join(out_dir,'rchie-WT-MUT')
    WT_hlx = csv_to_hlx(WT_csv,hybrid_energy=False, len_target=len_target, shift_len_query=shift_len_query)
    MUT_hlx = csv_to_hlx(MUT_csv,hybrid_energy=False, len_target=len_target, shift_len_query=shift_len_query)
    plot_rchie(WT_hlx, MUT_hlx,out_name, energy_large_low_mid_high=3,convert_png=convert_png,convert_svg=convert_svg)
    
    WT_hlx_hyb = csv_to_hlx(WT_csv,hybrid_energy=True, len_target=len_target, shift_len_query=shift_len_query)
    MUT_hlx_hyb = csv_to_hlx(MUT_csv,hybrid_energy=True, len_target=len_target, shift_len_query=shift_len_query)
    plot_rchie(WT_hlx_hyb, MUT_hlx_hyb,out_name+'.hybrid-energy',energy_large_low_mid_high=3,convert_png=convert_png,convert_svg=convert_svg)
    print("Plotted: ", out_name, out_name+'.hybrid-energy.pdf')
    return 

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
def get_seq_len(fa):
    return len(SeqIO.read(fa, "fasta"))
    
    
    
def run_intarna_plot_all(query_fasta_wt, target_fasta_wt, query_fasta_mut, target_fasta_mut, shift_len_query, out_dir='./'):
    
    len_target = get_seq_len(target_fasta_wt)
    
    csv_subopt_wt = os.path.join(out_dir, 'intarnaout_wt.csv')
    run_intaRNA_csv(query_fasta_wt, target_fasta_wt, "Q", "T", out_csv_file=csv_subopt_wt)

    csv_subopt_mut = os.path.join(out_dir, 'intarnaout_mut.csv')
    run_intaRNA_csv(query_fasta_mut, target_fasta_mut, "Q", "T", out_csv_file=csv_subopt_mut)


    plot_all(csv_subopt_wt, csv_subopt_mut, out_dir, len_target=len_target, shift_len_query=shift_len_query, convert_png=False)
    
# wt_mut_pdf = plot_all('../results/intarna/2019/WT-C100-C1000-5281772/predictions.sorted.csv', 
#                     '../results/intarna/2019/MUT-C100-C1000-7187161/predictions.sorted.csv', 
#                       './', 1000, convert_png=False, convert_svg=True)

# query_file, target_file = '../tmp/input_query.fa', '../tmp/input_target.fa'
# run_intaRNA_csv(query_file, query_file, "Q", "T", out_csv_file='../tmp/intarnaout.csv')
# plot_all('../tmp/intarnaout.csv', '../tmp/intarnaout.csv', '../tmp/', 100, convert_png=False)

def parse_SNP_tag(SNP_tag):
    #annot_locs += [int(SNP_tag[1:-1])]
    #annot_names += [SNP_tag
     
    matches =  re.match('(\D)(\d+)(\D)', SNP_tag)
    if not matches:
        raise RuntimeError("Invalid SNP tag: {}".format(SNP_tag)) 
    wild_char, loc, mut_char = matches.group(1), int(matches.group(2)), matches.group(3)
    return(wild_char, loc, mut_char)
                    
def rec_to_fa(rec, fasta_out):
    with open(fasta_out, 'w') as outfa:
        outfa.write(">{}\n{}\n".format(rec.id, rec.seq))
    

def run_all_SNP(fasta_wt, SNP_tag, query_context_len=20, out_dir='./'):
    
    # todo: check +-1 boundary cases of position and lenghts in this method
    
    wild_char, snp_loc, mut_char = parse_SNP_tag(SNP_tag)
    rec_wt = SeqIO.read(fasta_wt, "fasta")
    wild_seq = rec_wt.seq
    
    if len(wild_seq) < 200:
        print("NOTE!: Skipping long-range interaction predictions, sequence too short {} < 200".format(len(wild_seq)))
        return

    if wild_seq[snp_loc-1] != wild_char:
        print("WARNING!: SNP {} wild char expected: {}, but found non-matching:{} on wildtype sequences".format(SNP_tag, wild_char, wild_seq[snp_loc-1]))

    if query_context_len > snp_loc:
        print ("Warning: context len larger (={}) than possible context around SNP (={})\n context len shortened to {}".format(query_context_len, snp_loc, snp_loc-1))
        query_context_len = snp_loc-1

    if query_context_len > len(rec_wt)-snp_loc:
        print ("Warning: context len larger (={}) than possible context around SNP (={})\n context len shortened to {}".format(query_context_len, snp_loc, len(rec_wt)-snp_loc-1))
        query_context_len = len(rec_wt)-snp_loc -1 

    mut_seq = wild_seq[:snp_loc-1] + mut_char + wild_seq[snp_loc:]
    rec_mut = SeqRecord(mut_seq, id=rec_wt.id + '-MUTANT')
    fasta_mut = fasta_wt + '.mut.fa'
    rec_to_fa(rec_mut, fasta_mut)
    
    query_wt_seq = wild_seq[snp_loc-query_context_len:snp_loc+query_context_len]
    print("wild_seq is", wild_seq)
    print("snp_loc-query_context_len:snp_loc+query_context_len",snp_loc-query_context_len, snp_loc+query_context_len)
    print("query_wt_seq is", query_wt_seq)
    rec_query_wt = SeqRecord(query_wt_seq, id='{}-C{}'.format(rec_wt.id, query_context_len))
    fasta_query_wt = '{}.C{}.fa'.format(fasta_wt, query_context_len)
    rec_to_fa(rec_query_wt, fasta_query_wt)
    
    query_mut_seq = mut_seq[snp_loc-query_context_len:snp_loc+query_context_len]
    rec_query_mut = SeqRecord(query_mut_seq, id='{}-C{}-MUTANT'.format(rec_wt.id, query_context_len))
    fasta_query_mut = '{}.C{}.mut.fa'.format(fasta_wt, query_context_len)
    rec_to_fa(rec_query_mut, fasta_query_mut)
    
    shift_len = snp_loc - query_context_len 
    run_intarna_plot_all(fasta_query_wt, fasta_wt,
                         fasta_query_mut, fasta_mut,
                         shift_len_query=shift_len, out_dir=out_dir)
    
    
    
    

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='MutaRNA-long-range predict and plot long-range interactions of wildtype and mutant RNAs using IntaRNA'\
        '\nSample call: \"python MutaRNA-long-range.py --fasta-wildtype wt.fa --SNP-tag G3C --out-dir /tmp/ \"'
        )
    parser.add_argument('--fasta-wildtype', required=True, type=is_valid_file, help='Input sequence wildtype in fasta format')
    parser.add_argument('--SNP-tag',  required=True, type=str, help='SNP tag e.g. "C3G" for mutation at position 3 from C to G')
    parser.add_argument('--out-dir', default="./", type=is_valid_directory, help='output directory')
    parser.add_argument('--mutation-context-query',  default=20, type=int, help='Compute interamolecular interactions of the full sequence with this (short) context length around the mutation location. ')
    parser.add_argument('--rchie-bin-path', type=is_valid_file, help='Path to the rchie R script, e.g. ./r-chie/rchie.R')


    args = parser.parse_args()
    if args.rchie_bin_path:
        RCHIE_BIN = args.rchie_bin_path
    
    
    run_all_SNP(args.fasta_wildtype, args.SNP_tag, query_context_len=args.mutation_context_query, 
                out_dir=args.out_dir)
    
    
    