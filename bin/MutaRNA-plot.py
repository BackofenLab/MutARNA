#!/usr/bin/env python

import sys, os

script_path = os.path.dirname(__file__)
localdotplot_path = os.path.join(script_path, '../lib/local_dotplot/')
sys.path.append(localdotplot_path)
import local_dotplot_lib as ldp

snv_wrapper_path = os.path.join(script_path, '../lib/RNA_SNV_wrapper/')
sys.path.append(snv_wrapper_path)
import snv_wrapper


from matplotlib import pyplot as plt
import seaborn as sns

import numpy as np
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import itertools

from Bio.SeqRecord import SeqRecord
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import Alphabet
from subprocess import Popen, PIPE
from os.path import isfile

def is_valid_SNP(SNP_tag):
    SNP_tag = SNP_tag.strip()
    matches =  re.match('^(\D)(\d+)(\D)$', SNP_tag)
    
    if not matches:
        raise RuntimeError("Invalid SNP tag: \"{}\". A valid SNP tag example: G200C".format(SNP_tag))
    return SNP_tag

def is_valid_file(file_name):
    if os.path.isfile(file_name):
        return os.path.abspath(file_name)
    else:
        raise FileNotFoundError(os.path.abspath(file_name))
def is_valid_directory(dir_name):
    if os.path.isdir(dir_name):
        return os.path.abspath(dir_name)
    else:
        raise NotADirectoryError(os.path.abspath(dir_name))

def is_valid_sequence(s):
    rec = SeqRecord(Seq(s.upper().replace('T','U'), IUPAC.unambiguous_rna), id="RNA")
    if not Alphabet._verify_alphabet(rec.seq):
        raise RuntimeError("Invalid nucleotide sequence, unknown characters in input string {}".format(s))
    return rec
  
def makeSafeFilename(inputFilename):   
    try:
        safechars = string.letters + string.digits + " -_."
        return filter(lambda c: c in safechars, inputFilename)
    except:
        return ""  


def get_CD_record(full_rec, utr5_len, utr3_len):
    if utr3_len == 0:
        CD_seq = full_rec.seq[utr5_len:]
    else:
        CD_seq = full_rec.seq[utr5_len: -utr3_len]
    CD_rec = SeqRecord(CD_seq, id=full_rec.id+'-CDS', name=full_rec.name, description='utr5-0 utr3-0')
    return CD_rec


def write_dp_from_matrix(matrix, out_dp, template_dp,p_range=[-1,1]):
    from math import sqrt
    import re
    ureg = re.compile(r'^(\d+)\s+(\d+)\s+(\d+\.\d+)\s+[ul]box\s*')
    with open(template_dp) as template, open(out_dp, 'w') as out_handle:
        # Write first part before probs
        for line in template:
            if "ubox" in line or "lbox" in line:
                um = ureg.match(line)
                if um:
                    break
            out_handle.write(line)
        cx = matrix.tocoo()
        #Write alternative given probs
        for i,j,p in zip(cx.row, cx.col, cx.data):# used to have itertoolz.izip in Python2
            if abs(p) >= 1e-5 and p_range[0]<=p<=p_range[1]:
                out_handle.write("{} {} {:1.9f} ubox\n".format(i+1,j+1,sqrt(abs(p))))
        #Write remianing footer after probs
        for line in template:
            if "ubox" in line or "lbox" in line:
                um = ureg.match(line)
                if um:
                    continue
            out_handle.write(line)

def write_diff_dp(dp_wild, dp_mut, out_dp):
    '''Reads two dotplot matrices and write the absoloute difference of basepair probs into out_dp'''
    #     print "write_diff_dp: ", dp_wild, dp_mut, out_dp
    mut_matrix, mfe_dic = ldp.parse_dp_ps_sparse(dp_wild, sparse=True)
    wild_matrix, mfe_dic = ldp.parse_dp_ps_sparse(dp_mut, sparse=True)
    assert mut_matrix.shape == wild_matrix.shape
    diff_mat = (mut_matrix - wild_matrix)
    write_dp_from_matrix(diff_mat, out_dp=out_dp, template_dp=dp_wild)
    
    # dp_removed = '.'.join(
        # out_dp.split('.')[:-1]+['.removed.dp']+out_dp.split('.')[-1:])
    dp_removed = os.path.join(os.path.dirname(out_dp), os.path.basename(out_dp).replace('.ps','_weakened.dp'))
    write_dp_from_matrix(diff_mat, out_dp=dp_removed, 
                         template_dp=dp_wild,p_range=[0,1])
    # dp_introduced = '.'.join(
    #     out_dp.split('.')[:-1]+['.introduced.dp']+out_dp.split('.')[-1:])
    dp_introduced = os.path.join(os.path.dirname(out_dp), os.path.basename(out_dp).replace('.ps','_increased.dp'))

    write_dp_from_matrix(diff_mat, out_dp=dp_introduced, 
                         template_dp=dp_wild,p_range=[-1.0,0])
    return out_dp, dp_removed, dp_introduced 
    



def call_vienna_plfold(sequence, seq_name, do_localfold=False, local_W=200, local_L=150, global_L=1000, out_dir='./'):
    '''Runs Vienna RNAfold with partition function for all sequences inside input fasta file 
    # call_RNAfold_pf("ACCGGCUUAAAGG", "seq1")'''

    from subprocess import Popen, PIPE
    dp_file_name = os.path.join(out_dir, "{}_dp.ps".format(seq_name))
    unp_file_name = os.path.join(out_dir, "{}_lunp".format(seq_name))
    if isfile(dp_file_name): # Caution Race condition may occur 
        os.remove(dp_file_name)
    if isfile(unp_file_name): # Caution Race condition may occur 
        os.remove(unp_file_name)
    if not do_localfold:
        local_W, local_L = len(sequence), global_L
    
    import tempfile
    tmpfile = tempfile.NamedTemporaryFile(delete=False)

    # Open the file for writing.
    with open(tmpfile.name, 'w') as f:
        f.write('>{}\n{}\n'.format(seq_name, sequence)) 

    # RNAFOLD = 'RNAfold -p2 '
    RNAPLFOLD = 'RNAplfold -W {} -L {} -u 1 < {}'.format(local_W, local_L, tmpfile.name)  # -u 1 for unpaired probablitiy 
    assert len(sequence.split()) == 1
    cmd = "cd {}; ".format(out_dir)
    # cmd += ('echo ">%s\\n%s\\n" | '%(seq_name, sequence))  
    cmd += RNAPLFOLD 
    # cmd += '; ls'
    print(cmd)
    p = Popen(cmd, stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    print (out)
    print (err)
    has_error = False
    if err:
        if b"warn" in err.lower():
            print ("Warning in calling call_RNAfold:\n {} {}\n".format(out, err))
        else:
            raise RuntimeError("Error in calling call_RNAfold: {} {}\n".format(out, err))
    if not isfile(dp_file_name):
        raise RuntimeError("Error: Expected dp file: {} is not created!".format(dp_file_name))
    if not isfile(unp_file_name):
        raise RuntimeError("Error: Expected lunp file: {} is not created!".format(unp_file_name))
    
    return dp_file_name, unp_file_name

    
def run_dot2circ(dp_file, prefix, out_dir=""):
    file_name_string = "".join(x for x in dp_file if x not in ['|','<', '>'])
    #print ("dp_file:", dp_file, file_name_string)
    
    cmd = 'cd "{}/dot2circ/"; python dot2circ.py --outputdir \'{}\' --prefix \'{}\' --dp-file \'{}\' --title \'{}\' '.format(
        script_path, out_dir, prefix+'-circos', file_name_string, prefix)
    print (cmd)
    p = Popen(cmd, stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    print (out.decode('ascii'))
#     print (err.decode('ascii'))
    has_error = False
    if err:
        raise RuntimeError("Error in calling dot2circ.py: {} {}\n".format(out.decode('utf-8'), err.decode('utf-8')))

    
def create_circos_annotation(CDS_len, utr5_len, utr3_len, snp_locs, snp_names, CDS_label='RNA'):
    '''
    Genes formatted Example:
    seq 193 759 CDS fill_color=green,r0=1.01r,r1=1.01r+20p
    seq 0 192  5UTR fill_color=yellow,r0=1.01r,r1=1.01r+20p
    seq 760 5765  3UTR fill_color=yellow,r0=1.01r,r1=1.01r+20p
    seq 227 228 SN-86 fill_color=vdred,r0=1.03r,r1=1.03r+20p
    '''
    print ("snp_locs:", snp_locs)
    start = 0
    formatted_str = ""
    if utr5_len is not None:
        formatted_str += 'seq {} {} 5` fill_color=yellow,r0=1.01r,r1=1.01r+20p\n'.format(start, start+utr5_len)
        start += utr5_len + 1 # Tocheck: maybe plus one not needed?
    if CDS_len is not None:
        
        formatted_str += 'seq {} {} {} fill_color=green,r0=1.01r,r1=1.01r+20p\n'.format(start, start+CDS_len,  CDS_label)
        start += CDS_len + 1 # Tocheck: maybe plus one not needed?
    if utr3_len is not None:
        formatted_str += 'seq {} {} 3` fill_color=blue,r0=1.01r,r1=1.01r+20p\n'.format(start, start+utr3_len)
        start += utr3_len + 1 # Tocheck: maybe plus one not needed?
    if snp_locs is not None:
        for i in range(len(snp_locs)):
            formatted_str += 'seq {} {} c.-{} fill_color=vdred,r0=1.03r,r1=1.03r+20p\n'.format(snp_locs[i], snp_locs[i]+1, snp_names[i])
    # print (formatted_str)
    genes_file =  '{}/dot2circ/data/genes.formatted.txt'.format(script_path)
    with open (genes_file, 'w') as genes_out:
        genes_out.write(formatted_str)

def plot_circos_seq_annotate(rec, annotate_indices, annotate_names, local_fold=False, 
                             plotted_seq_lenght=None,utr5_l=0, utr3_l = 0,color='r',dp_full = None,suffix='' ):
 
    seq = rec.seq

    if len(annotate_indices) != len(annotate_names):
        raise RuntimeError("Mismatch indic/names: {}, {} }".format(len(annotate_indices), len(annotate_names))) 
    
    ID = rec.id
    if dp_full is None:
        dp_full, unp_full = call_vienna_plfold(rec.seq, ID, local_fold)
    create_circos_annotation(len(rec)-utr5_l-utr3_l, utr5_l, utr3_l, annotate_indices, annotate_names)

    run_dot2circ(dp_full, ID+suffix)
        
    #dp_diff = dp_mut.replace('.ps', '_diff.ps')
    #write_diff_dp(dp_full, dp_mut, dp_diff)
    #run_dot2circ(dp_diff, rec_mut.id+suffix+'-diff')
    

from cycler import cycler



def get_unpaired_probs(unp_file):
    # Read Vienna RNAplfold unpaired prob file (-u 1) into dict
    with open(unp_file) as unp_in:
        line1 = unp_in.readline()
        if "#unpaired probabilities" not in line1:
            raise IOError('Unexpected header for lunp file: {}'.format(line1))
        line2 = unp_in.readline()
        if "#i$\tl=1" not in line2:
            raise IOError('Unexpected second header for lunp file: {}'.format(line2))
        up_dic = dict()
        for line in unp_in:
            splits = line.split()
            assert len(splits) >= 2
            pos = int(splits[0])
            up_prob = float(splits[1])
            assert pos >=1
            assert up_prob >= 0 and up_prob <= 1
            assert pos not in up_dic
            up_dic[pos] = up_prob
            
    #print('Parsed ', unp_file)
    return up_dic


def plot_up_dict(up_dic, plot_lims=None, title='XX', fig=None, diff=False,tidy=False):
    if plot_lims is None:
        x, y = up_dic.keys(), up_dic.values()
    else:
        x, y = up_dic.keys()[plot_lims[0]:plot_lims[1]+1], up_dic.values()[plot_lims[0]:plot_lims[1]+1]
    if fig is None:
        fig = plt.figure(figsize=(9, 3))
    x, y = list(x), list(y)
    # print(list(zip(x,y)))
    if fig.get_axes():
        ax = fig.get_axes()[0]
    else:
        ax = fig.add_subplot(111) 
    ax.plot(x, y, label=title, alpha=0.8)
    if not tidy:
        ax.legend(loc='lower left', framealpha=0.2)#, bbox_to_anchor=(0.0, 1.1))



    if len(x) < 101:
        ticks_label_step = 10
        ticks_step = 10
    else:
        ticks_label_step = 50
        ticks_step = 10

    minor_ticks = np.arange(min(x), max(x), 1)                                               
    major_ticks = np.arange(min(x)-min(x)%ticks_step, max(x), ticks_step)                                               

    ax.set_xticks(minor_ticks, minor=True)                                           
    ax.set_xticks(major_ticks) 
                

    ax.set_xlabel('Position')
    if diff:
        ax.set_yticks([-1, 1], minor=True)                                           
        ax.set_ylim([-1.05,1.05])
        #ax.set_ylabel('Accessibility(WT) - Accessibility(MUT)')
    else:
        ax.set_yticks([0,1], minor=True)                                           
        ax.set_ylim([-0.05,1.05])
    ax.set_ylabel('Accessibility')
#     ax.grid(which='both')                                                            

    # or if you want different settings for the grids:                               
    ax.grid(which='minor', alpha=0.5)
    ax.axhline(0)
#     ax.axhline(0, linestyle='--', color='k', alpha=0.5) # horizontal lines
#     ax.axhline(1, linestyle='--', color='k', alpha=0.5) # horizontal lines
    
    ax.set_xlim([min(x)-1, max(x)+1])
    # if not tidy:
        # ax.set_title(title)
    
    sort_legends = True
    if sort_legends:
        handles, labels = ax.get_legend_handles_labels()
        # sort both labels and handles by labels
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
        ax.legend(handles, labels)


    
    if ticks_label_step != ticks_step:
        labels = [item.get_text() for item in ax.get_xticklabels()]
        labels_locs = ax.get_xticks()
        pruned_labels = [str(loc)  if (loc%ticks_label_step)==0 else '' for loc, lab in zip(labels_locs, labels)]
        ax.set_xticklabels(pruned_labels)



def heatmap_up_dict(up_dic, plot_lims=None, title='XX', ax=None,fig=None, diff=False,
                    ticklabel=True,legend=True,ticks_step=10):
    if plot_lims is None:
        x, y = list(up_dic.keys()), list(up_dic.values())
    else:
        x, y = list(up_dic.keys())[plot_lims[0]:plot_lims[1]+1], list(up_dic.values())[plot_lims[0]:plot_lims[1]+1]
    print(len(x), len(y))
    if(len(x)>200):
        print("Skipping heatmap, sequence longer than 200 limit. < {} ".format(len(x)))
        return
    
    if ax is None:
        ax = fig.add_subplot(121)
    ax.yaxis.set_label_position("right")

    # Now adding the colorbar
#     pos1 = ax.get_position() # get the original position 
#     pos2 = [pos1.x0 + 0.3, pos1.y0 + 0.3,  pos1.width / 2.0, pos1.height / 2.0]
    if legend:
        cax = fig.add_axes([1.0, 0.5, 0.1, .1])
    else:
        cax = None

#     ax.plot(x, y, label=title)
#     sns.pointplot(x,y,ax=ax)
    ax.set_xticks(np.arange(min(x)-min(x)%ticks_step, max(x), ticks_step))                                               
    ax.set_yticks(np.arange(min(y)-min(y)%ticks_step, max(y), ticks_step))                                              

     

    sns.heatmap(np.reshape(np.array(y),(len(y),1)), cmap=sns.cubehelix_palette(as_cmap=True), ax=ax, 
                cbar=legend,cbar_ax=cax,vmin=0.0, vmax=1.0,
                # yticklabels=x, 
               )
#     cbar = ax.collections[0].colorbar
#     cbar.set_ticks([0., .2, .4, .6, .8, 1.0])
    
#     cbar.set_ticklabels(['low', '20%', '75%', '100%'])
    if legend:
        ax.legend(loc='upper left')#, bbox_to_anchor=(0.0, 1.1))

    minor_ticks = np.arange(min(x), max(x), 1)                                               
    major_ticks = np.arange(min(x)-min(x)%10, max(x), 10)                                               
    yticks = x
    keptticks = yticks[::10]
    yticks = ['' for y in yticks]
    yticks[::10] = keptticks
#     ax.set_yticks(minor_ticks,                   minor=True)
    from matplotlib.ticker import AutoMinorLocator,MultipleLocator,FormatStrFormatter,FuncFormatter,IndexLocator,StrMethodFormatter
    minorLocator = IndexLocator(1,offset=0.5)
    majorLocator = IndexLocator(10,offset=0.5)
    def incer(x, pos):
        'The two args are the value and tick position'
        return '%d' % (x+1)

    majorFormatter = FuncFormatter(incer)
#     majorFormatter = StrMethodFormatter('{x}',use_offset=False                                       )

#     ax.set_yticklabels(yticks, rotation=0,)

    ax.yaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_major_locator(majorLocator)
    ax.tick_params(which='major', length=6)
    ax.tick_params(which='minor', length=3)
    if ticklabel:
        ax.yaxis.set_major_formatter(majorFormatter)
#         ax.set_yticklabels(yticks, rotation=0,)

#         ax.tick_params(which='minor', length=6)

        ax.yaxis.set_label_position("right")
    else:
        ax.set_yticklabels([])
#     ax.set_yticks(major_ticks) 
#     ax.set_ylabel('Position')
#     ax.set_yticks([-1, 1], minor=False, )                                           

    if diff:
        ax.set_xticks([-1, 1], minor=True)                                           
        ax.set_xlim([-1.05,1.05])
        ax.set_xlabel('P_unpaired(wild) - P_unpaired(mut)')
    else:
        ax.set_xticks([], minor=False,)                                           
#         ax.set_xlim([-0.05,1.05])
#         ax.set_xlabel('Accessibility')
#     ax.grid(which='both')                                                            

    # or if you want differnet settings for the grids:                               
#     ax.grid(which='minor', alpha=0.5) 
#     ax.axhline(0, linestyle='--', color='k', alpha=0.5) # horizontal lines
#     ax.axhline(1, linestyle='--', color='k', alpha=0.5) # horizontal lines
    
#     ax.set_ylim([min(x)-1, max(x)+1])
    ax.set_title(title,rotation=90,va='bottom')


def plot_unpaired_probs(up_file_pairs, plot_heatmap=False,rang=None, out_dir='./',ECGs_together=True, dynamic_width_ecg=True):
    plt.rc('axes', #prop_cycle=(cycler('color', ['k',(0.8901960784313725, 0.10196078431372549, 0.10980392156862745)]
            prop_cycle=(cycler('color', sns.color_palette('colorblind', 3))))
            
            #           sns.color_palette()
            #           [sns.color_palette("Paired")[1], sns.color_palette("Paired")[-1]]
            #           #                                   sns.color_palette("RdBu_r", 2)
            #           ['k', 'b', 'r', 'g']
            #          ) 
            #    +                            cycler('linestyle', ['-', '--', ':', '-.']) 
            #   ))

   
    for up_file_wild, up_file_mut in up_file_pairs:
        # print (up_file_wild, up_file_mut)
        d_wild = get_unpaired_probs(up_file_wild)
        d_mut = get_unpaired_probs(up_file_mut)
        dict_diff = {f:(d_wild[f] - d_mut[f])  for f in d_wild}
        seq_len = len(d_wild)
        if dynamic_width_ecg:
            fig_width = 5 + 2*int(seq_len/100)
        else:
            fig_width = 12

        if plot_heatmap:
            title_key = 'heatband'
            fig, ax = plt.subplots(1,figsize=(0.5, 11))
            #heatmap_up_dict(dict_diff, rang, title=os.path.basename(up_file_wild).replace('-WT-','-').replace('_lunp','-DIFF'), ax=ax, fig=fig,diff=True)
        else:
            title_key = 'ECG'

            if ECGs_together:
                fig = plt.figure(figsize=(fig_width, 3))
            else:
                fig = plt.figure(figsize=(fig_width, 2))
            plot_up_dict(dict_diff, rang, title='Accessibility(wt) - Accessibility(mut)' #os.path.basename(up_file_wild).replace('-WT-','-').replace('_lunp','-DIFF')
                         , fig=fig, diff=True,tidy=True)
        fig.savefig(os.path.join(out_dir, os.path.basename(up_file_wild)+'-diff-{}.png'.format(title_key)), bbox_inches='tight', pad_inches=0.2
                   )
        fig.savefig(os.path.join(out_dir, os.path.basename(up_file_wild)+'-diff-{}.svg'.format(title_key)), bbox_inches='tight', pad_inches=0.2
                   )
        
        if plot_heatmap:
            fig = plt.figure(figsize=(2.5, 11),tight_layout={'w_pad':2})

        elif not ECGs_together:
            fig = plt.figure(figsize=(fig_width, 1))
        labeldic = {1:True, 0:False}
        titledic = {0:'mut', 1:'wt'}
        
        for iup, up_file in enumerate([up_file_mut, up_file_wild]):
            if plot_heatmap:
                ax = fig.add_subplot(131+iup*2) # Skip one ax in dirty way for cbar
                heatmap_up_dict(get_unpaired_probs(up_file), rang, title=os.path.basename(up_file).replace('-MUT-','-').replace('-WT-','-').replace('_lunp',''), 
                                fig=fig,ax=ax,ticklabel=labeldic[iup])            
            else:
                plot_up_dict(get_unpaired_probs(up_file), rang, title='Accessibility({})'.format(titledic[iup]),#os.path.basename(up_file).replace('-MUT-','-').replace('-WT-','-').replace('_lunp',''), 
                             fig=fig,tidy=False,diff=ECGs_together
#                              ax=ax, ticklabel=labeldic[iup]
                            )            
            
#         fig.tight_layout(pad=0.1)

        fig.savefig(os.path.join(out_dir, os.path.basename(up_file_wild).replace('WT','WTMUT')+'-{}.png'.format(title_key)), bbox_inches='tight', #pad_inches=0.5,
                    dpi=600)
        fig.savefig(os.path.join(out_dir, os.path.basename(up_file_wild).replace('WT','WTMUT')+'-{}.svg'.format(title_key)), bbox_inches='tight', #pad_inches=0.5,
                    )

def plot_circos_seq_SNP(rec_wild, SNP_tag, rec_mut, do_local=True,do_global=False, plotted_seq_lenght=None,
dotplot=True,ECGplot=True,suffix='',annot_locs=[], annot_names=[],local_global_out_dir='./', local_L=150, local_W=200, global_L=1000):
    
    ID = '_'.join(rec_wild.id.split('|')[:2]) # +'_'+SNP_tag#"".join(x for x in rec_wild.id if x not in ['|','<', '>'])
    utr5_l, utr3_l = 0, 0 
    #     print "rec_wild: " , rec_wild

    local_fold_runs = []
    if do_local:
        ldir = os.path.join(local_global_out_dir, 'local/')
        os.makedirs(ldir, exist_ok=True) 
        local_fold_runs   += [(True,ldir)] 
    if do_global:
        gdir = os.path.join(local_global_out_dir, 'global/')
        os.makedirs(gdir, exist_ok=True) 
        local_fold_runs   += [(False,gdir)] 

    for (local_fold, out_dir) in local_fold_runs:
        dp_wild, unp_wild = call_vienna_plfold(rec_wild.seq, ID, local_fold, local_L=local_L, local_W=local_W, global_L=global_L, out_dir=out_dir)
        create_circos_annotation(len(rec_wild), utr5_l, utr3_l, annot_locs, annot_names)
        run_dot2circ(dp_wild, ID+'-WILD'+suffix, out_dir=out_dir)
        

        dp_mut, unp_mut = call_vienna_plfold(rec_mut.seq, rec_mut.id, local_fold, local_L=local_L, local_W=local_W, global_L=global_L,out_dir=out_dir)

        
        if len(SNP_tag) > 0 :
            matches =  re.match('(\D)(\d+)(\D)', SNP_tag)
            if not matches:
                raise RuntimeError("No matches founs for tag:{}".format(SNP_tag)) 
            wild_char, loc, mut_char = matches.group(1), int(matches.group(2)), matches.group(3)

            annot_locs += [loc]
            annot_names += [SNP_tag]

        create_circos_annotation(len(rec_mut), utr5_l, utr3_l, annot_locs, annot_names)
        run_dot2circ(dp_mut, rec_mut.id+suffix, out_dir=out_dir)
        dp_diff = dp_mut.replace('.ps', '_diff.ps')
        dpabs, dpremove, dpintroduce = write_diff_dp(dp_wild, dp_mut, dp_diff)

        #run_dot2circ(dp_diff, rec_mut.id+suffix+'-diff', out_dir=out_dir)
        run_dot2circ(dpremove, rec_mut.id+suffix+'-weakened', out_dir=out_dir)
        run_dot2circ(dpintroduce, rec_mut.id+suffix+'-increased', out_dir=out_dir)
        
        if dotplot is True:
            #ldp.plot_heat_maps(None, ldp.parse_dp_ps(dp_wild), filename=ID+'-WILD', title_suffix=ID+'\n'+r'$P({\rm WT})$', what='basepairs',inverse=True, out_dir=out_dir)#, gene_loc=[2,10])
            #ldp.plot_heat_maps(None, ldp.parse_dp_ps(dp_mut), filename=ID+'-MUTANT', title_suffix=ID+'\n'+r'$P({\rm mutant})$''\n'+r'$P({\rm wt})$''-MUTANT', what='basepairs',inverse=True, out_dir=out_dir)

            ldp.plot_heat_maps(None, ldp.parse_dp_ps(dp_wild)-ldp.parse_dp_ps(dp_mut), colormap='seismic', vmin=-1.0, vmax=1.0,
                                filename=ID+'-DIFF',title_suffix=ID+'\n'+r'$\Delta = P({\rm WT})-P({\rm mutant})$', what='basepairs',inverse=True, out_dir=out_dir)
            #ldp.plot_heat_maps(None, ldp.parse_dp_ps(dp_diff), filename=ID+'-ABSDIFF',title_suffix=ID+'-ABSDIFF', what='basepairs',inverse=True, out_dir=out_dir)
            #ldp.plot_heat_maps(None, ldp.parse_dp_ps(dpremove), filename=ID+'-REMOVED', title_suffix=ID+'-REMOVED', what='basepairs',inverse=True, out_dir=out_dir)
            #ldp.plot_heat_maps(None, ldp.parse_dp_ps(dpintroduce), filename=ID+'-INTRODUCED', title_suffix=ID+'-INTRODUCED', what='basepairs',inverse=True, out_dir=out_dir)
            
            ldp.plot_heat_maps(None, ldp.parse_dp_ps(dp_wild)+ldp.parse_dp_ps(dp_mut).transpose(), filename=ID+'-WT-MUT', what='basepairs',
                    inverse=True, interactive=False, gene_loc=None,title_suffix=ID+'-'+SNP_tag+'\n'r'$P({\rm WT})$, $P({\rm mutant})$', out_dir=out_dir, upper_triangle_txt='wt',lower_triangle_txt='mut')
            
            ldp.plot_heat_maps(None, ldp.parse_dp_ps(dpremove)+ldp.parse_dp_ps(dpintroduce).transpose(), filename=ID+'-REMOVED-INTRODUCED', what='basepairs',
                    inverse=True, interactive=False, gene_loc=None,title_suffix=ID+'\n'+r'$|\Delta| = |P({\rm WT})-P({\rm mutant})|$', out_dir=out_dir, upper_triangle_txt='weakened\n' + r'    $\Delta>0$',lower_triangle_txt='increased\n' + r'    $\Delta<0$')

        if ECGplot is True:
    #         plot_up_dict(u, None, title=ID, fig=myfig,tidy=True)
            plot_unpaired_probs([(unp_wild, unp_mut)], plot_heatmap=False, out_dir=out_dir)
            #plot_unpaired_probs([(unp_wild, unp_mut)], plot_heatmap=True, out_dir=out_dir)

def get_mutation_rec(wild_rec, SNP_tag):
    wild_seq = wild_rec.seq
    matches =  re.match('(\D)(\d+)(\D)', SNP_tag)
    if not matches:
        raise RuntimeError("No matches founs for tag:{}".format(SNP_tag)) 
    wild_char, loc, mut_char = matches.group(1), int(matches.group(2)), matches.group(3)
    if len(wild_seq) < loc:
        raise RuntimeError("SNP loc outside sequence len:{}".format(SNP_tag)) 

    if (wild_seq[loc-1].upper() != wild_char.upper()):
        print("WARNING!: SNP {} wild char expected: {}, but found non-matching:{} on wildtype sequences".format(SNP_tag, wild_char, wild_seq[loc-1]))

    mut_seq = wild_seq[:loc-1] + mut_char + wild_seq[loc:]
    
    #print mut_seq
    rec_mut = SeqRecord(mut_seq, id=wild_rec.id + '-MUTANT')
    return rec_mut

def filter_SNV_columns(df, clean_columns=None, select_not_filter=True):
    if (select_not_filter is True) and (clean_columns is not None):
        print(clean_columns)
        clean_columns += ['tool']
        return df[clean_columns].copy()
    if clean_columns is None:
        clean_columns = ['tool','SNP', 'd', 'd_max', 'interval', 'interval.1',
        'p-value', 'p-value.1', 'r_min', 'rnasnp_params', 'w', 
        'MFE(wt)', 'MFE(mu)', 'dMFE', 'H(wt||mu)']
    clean_columns += ['tool', 'rnasnp_params']
    return df.loc[:, df.columns.isin(clean_columns)].copy()


def get_SNV_scores(fasta_wt, SNP_tag, out_dir='./'):

    df_remuRNA = snv_wrapper.run_remuRNA(fasta_wt, [SNP_tag], window=None)
    df_remuRNA['tool'] = 'remuRNA'
    df_remuRNA = filter_SNV_columns(df_remuRNA, ['SNP','H(wt||mu)', 'MFE(mu)', 'MFE(wt)', 'dMFE'])

    df_RNAsnp1 = snv_wrapper.run_RNAsnp(fasta_wt, [SNP_tag], window=None, plfold_W=None, plfold_L=None, mode=1)
    df_RNAsnp1['tool'] = 'RNAsnp'
    df_RNAsnp1 = filter_SNV_columns(df_RNAsnp1, ['SNP','interval', 'd_max', 'p-value']).rename(columns={'d_max':'distance', 'rnasnp_params':'params'})

    df_RNAsnp2 = snv_wrapper.run_RNAsnp(fasta_wt, [SNP_tag], window=None, plfold_W=None, plfold_L=None, mode=2)
    df_RNAsnp2['tool'] = 'RNAsnp'
    df_RNAsnp2 = filter_SNV_columns(df_RNAsnp2, ['SNP', 'interval', 'd',  'p-value']).rename(columns={'d':'distance', 'rnasnp_params':'params'})


    df_RNAsnp12 = pd.concat([df_RNAsnp1, df_RNAsnp2], sort=True)
    
  
    
    csv_remuRNA = os.path.join(out_dir, 'remuRNA.csv')
    csv_RNAsnp1 = os.path.join(out_dir, 'RNAsnp_mode1.csv')
    csv_RNAsnp2 = os.path.join(out_dir, 'RNAsnp_mode2.csv')
    csv_RNAsnp12 = os.path.join(out_dir, 'RNAsnp.csv')


    df_remuRNA.to_csv(csv_remuRNA, index=False)
    df_RNAsnp1.sort_index(axis=1).to_csv(csv_RNAsnp1, index=False)
    df_RNAsnp2.sort_index(axis=1).to_csv(csv_RNAsnp2, index=False)
    df_RNAsnp12.sort_index(axis=1).to_csv(csv_RNAsnp12, index=False)
    print('SNP scores were saved to:', csv_remuRNA, csv_RNAsnp12)
    return (csv_remuRNA, csv_RNAsnp12)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='MutaRNA-plot predict and plot local and global base-pair probabilities of wildtype and mutant RNAs'\
        '\nSample call: \"python bin/MutaRNA-plot.py --fasta-wildtype data/sample0.fa --SNP-tag G3C --out-dir tmp --no-global-fold\"'
        )
    
    parser.add_argument('--fasta-wildtype', required=True, type=is_valid_file, help='Input sequence wildtype in fasta format')
    #parser.add_argument('--sequence-wild', required=True, type=is_valid_sequence, help='Input sequence string wildtype')
    #parser.add_argument('--sequence-mutant',type=is_valid_sequence, help='Input sequence string mutant (support disabled)')
    parser.add_argument('--SNP-tag',  required=True, type=is_valid_SNP, help='SNP tag e.g. "C3G" for mutation at position 3 from C to G')
    parser.add_argument('--out-dir', default="./", type=is_valid_directory, help='path the output directory. The directory must already exist.')
    parser.add_argument('--no-global-fold', action='store_true', help='Do not run (semi-)global fold (semi: max-window 1000nt)')
    parser.add_argument('--no-local-fold', action='store_true', help='Do not run local fold')
    parser.add_argument('--local-W',  default=200, type=int, help='Window length for local fold')
    parser.add_argument('--local-L',  default=150, type=int, help='Max base-pair interaction span for local fold')
    parser.add_argument('--global-maxL',  default=1000, type=int, help='Maximum interaction span of global length.')
    parser.add_argument('--no-SNP-score', action='store_true', help='Do not run SNP structure abberation scores with RNAsnp and remuRNA')
    parser.add_argument('--enable-long-range', action='store_true', help='predict and plot long-range interactions of wildtype and mutant RNAs using IntaRNA')
    parser.add_argument('--enable-global-fold', action='store_true', help='enable global fold')



# Save to file in the current working directory

    args = parser.parse_args()
    
    #rec_wild = args.sequence_wild
    rec_wild = SeqIO.read(args.fasta_wildtype, 'fasta')
    rec_wild.id = "RNA"

    args.sequence_mutant = None # Disable sequence option 

    if args.sequence_mutant is None and args.SNP_tag is None:
        raise RuntimeError("Exactly one of these options must be passed (--sequence-mutant, --SNP-tag) but none is provided.")
    
    if args.sequence_mutant is not None and args.SNP_tag is not None:
        raise RuntimeError("Exactly one of these options must be passed (--sequence-mutant, --SNP-tag) but both are provided.")
    
    if args.sequence_mutant is None:
        rec_mutant = get_mutation_rec(rec_wild, args.SNP_tag)
        SNP_tag = args.SNP_tag
    else:
        rec_mutant = args.sequence_mutant
        rec_mutant.id = rec_wild.id + '-MUTANT'
        SNP_tag = ""

    if args.local_L > len(rec_wild):
        print ("Note: global and local outputs would be the same, since sequence length is shorter than bp-interaction lengnth. ")
        #raise RuntimeError ("Wildtype and mutant sequences have unequal lengths. wild:{} != mutant:{}".format(len(rec_mutant), len(args.sequence_wild)))

   
    plot_circos_seq_SNP(rec_wild, SNP_tag, rec_mut=rec_mutant, do_local=not args.no_local_fold, do_global=(not args.no_global_fold) and (args.enable_global_fold), 
    local_global_out_dir=args.out_dir, local_L=args.local_L, local_W=args.local_W, global_L=args.global_maxL)
    
    if not args.no_SNP_score:
        get_SNV_scores(args.fasta_wildtype, SNP_tag, out_dir=args.out_dir)

    if args.enable_long_range is True:
        lrdir = os.path.join(args.out_dir, 'long-range/')
        os.makedirs(lrdir, exist_ok=True)     
        import MutaRNA_long_range as mutLR

        mutLR.run_all_SNP(args.fasta_wildtype, args.SNP_tag,
                out_dir=lrdir)


