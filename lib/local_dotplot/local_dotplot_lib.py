import RNA
import numpy as np 
from altschulEriksonDinuclShuffle import dinuclShuffle
import random
import os
import re

from scipy.sparse import csr_matrix
import warnings
from scipy.sparse import SparseEfficiencyWarning
warnings.simplefilter('ignore', SparseEfficiencyWarning)

kT = (37+273.15)*1.98717/1000.0  # /* kT in kcal/mol */

VIENNA_BIN_PATH = '' #'/home/milad/1software/bin/'
RNAFOLD = 'RNAfold -p --noPS  '
RNAPLFOLD = 'RNAplfold '
QUAKE_PARAM_FILE = '/home/milad/1workspace/mmfold/src/misc/rna_turner2004_ML_up_penalty.par '
ANDERO_PARAM_FILE = '/home/milad/1workspace/mmfold/src/misc/rna_andronescu2007.par '


def compute_part_func(infile_fa, seq_names, outdir_path="./", use_plfold=False, which_params='turner',  dangles=2,
                      use_cache=False):
    '''Runs Vienna RNAfold/RNAplfold with partition function for all sequences inside input fasta file
    If use_cache, it does nothing if If the ps file with same paramaters exists '''
    from subprocess import Popen, PIPE
    #     print "compute_part_func(", infile_fa, seq_names
    if use_plfold:
        out_dir = outdir_path + "/" + RNAPLFOLD.replace(' ', '')
    else:
        out_dir = outdir_path + "/" + RNAFOLD.replace(' ', '')
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    if not os.path.isfile(infile_fa):
        raise IOError("Fastafile not found: {}".format(infile_fa))

    all_in_cache = all([os.path.isfile(os.path.join(out_dir, sname+'_dp.ps')) for sname in seq_names])
    if all_in_cache and use_cache:
        raise NotImplementedError("Sequence names for caching are not correctly set")
        return out_dir

    with open(infile_fa) as in_rna:
        arg_param = ""
        if which_params == 'quake':
            arg_param += " -P %s " % QUAKE_PARAM_FILE
#             dangles = 0
        elif which_params.startswith('andero'):
            assert dangles == 2
            arg_param += " -P %s " % ANDERO_PARAM_FILE
        elif which_params != 'turner':
            assert dangles == 2
            raise RuntimeError("Unknown parameter option {}".format(which_params))

        assert dangles >= 0 and dangles <= 2
        arg_param += " --dangles {} ".format(dangles)

        if use_plfold:
            p = Popen(('cd %s;' %out_dir) + VIENNA_BIN_PATH + RNAPLFOLD + arg_param, stdin=in_rna, shell=True, stdout=PIPE, stderr=PIPE)
        else:
            p = Popen(('cd %s;' %out_dir) + VIENNA_BIN_PATH +  RNAFOLD + arg_param, stdin=in_rna, shell=True, stdout=PIPE, stderr=PIPE)

        out, err = p.communicate()
        if err:
            print ("Error in calling RNAfold for ", infile_fa)
            print (out)
            print (err)

            # With long sequences RNAfold prints scalign factor to stderr
            if (not use_plfold and not ("scaling factor" in err or "free energy" in err)):
                raise RuntimeError

    return out_dir


def compute_mfe_probability(in_seq):
    '''Runs Vienna RNAfold/RNAplfold with partition function for all sequences inside input fasta file
    If use_cache, it does nothing if If the ps file with same paramaters exists '''

    from subprocess import Popen, PIPE

    p = Popen(('echo "%s" | ' % in_seq) + VIENNA_BIN_PATH + RNAFOLD, stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)

    out, err = p.communicate()
    if err:
        print ("Error in calling RNAfold for ", in_seq)
        print (out)
        print (err)
        raise RuntimeError
    lines = out.split('\n')
    #     print out
    mfe_line = lines[1]  # mfe line
    assert len(mfe_line.split()) >= 2
    mfe_energy = float(mfe_line.split()[-1].replace(')', '').replace('(', ''))

    pf_line = lines[2]  # pf line
    assert len(pf_line.split()) >= 2
    pf_energy = float(pf_line.split()[-1].replace(']', '').replace('[', ''))

    import re
    my_re = re.compile(r'\s*frequency of mfe structure in ensemble\s+([\d\.\-e]+);\sensemble diversity\s+([\d\.\-]+)')
    freq_line = lines[4]
    match = my_re.match(freq_line)
    if match is None or match.groups() is None:
        print ("Error unexpected frequency line format, found:\n  {}\n  {}\n".format(in_seq, freq_line))
        raise RuntimeError
    mfe_prob = (float)(match.groups()[0])
    diversity = (float)(match.groups()[1])
    return mfe_energy, pf_energy, mfe_prob, diversity


def getBPPM(sequence, structure="", bppm_cutoff = 0):
    import numpy as np
    seq_len = len(sequence)

    bppm = np.zeros((seq_len, seq_len))

    RNA.pf_fold(sequence, structure)
    for l_pos in xrange(0, seq_len):
        for r_pos in xrange(l_pos, seq_len+1):
            if l_pos<r_pos:
                bpp = RNA.get_pr(l_pos+1, r_pos+1)
                if bpp > bppm_cutoff:
                    bppm[l_pos, r_pos] = bpp
    RNA.free_pf_arrays()
    # print bppm
    return bppm


def get_mfe_probs(rna_context_seq):

    structure = ""
    structure, part_funct_whole = RNA.pf_fold(rna_context_seq, structure)
    RNA.free_arrays()
    structure, mfe_whole = RNA.fold(rna_context_seq, structure)
    RNA.free_arrays()

    from math import exp
    # print mfe_whole, part_funct_whole, exp(mfe_whole/(-kT))/exp(part_funct_whole/(-kT))
    # full_seq_outer = rna_seq[:rna_split_pos] +  selection + rna_seq[rna_split_pos:]
    # for context_len in range (0:len(context), 10):

    count = 0
    mfe_probs = np.zeros((len(rna_context_seq), len(rna_context_seq)))

    for l_pos in range(0, len(rna_context_seq)):
        for r_pos in range(l_pos+3, len(rna_context_seq)):
            sub_seq = rna_context_seq[l_pos:r_pos+1]
            bp_prob = RNA.get_pr(l_pos+1, r_pos+1)
            if bp_prob < 1e-3:  # 0:
                continue
    #       print l_pos, r_pos, bp_prob
            mfe_subseq, pf_subseq, mfe_prob_subseq, diversity_subseq = compute_mfe_probability(sub_seq)
            mfe_probs[l_pos, r_pos] = mfe_prob_subseq
            count += 1
    print (count)
    return mfe_probs


# ###########################################################################################
# Plotting visualiztion

import matplotlib.pyplot as plt
import numpy as np
import matplotlib

# import matplotlib as mpl
# mpl.rc("savefig", dpi=600)


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def my_heatmap_imshow(mat, fig, ax, threshold=1e-3, inverse=True, interactive=False):
    im = ax.imshow(mat, extent=(10, 20, 10, 20), origin='lower', zorder=1, interpolation='none')

#     plugins.connect(fig, plugins.MousePosition(fontsize=14))


def my_heatmapMatshowSparse(mat, fig, ax, threshold=1e-3, inverse=True, interactive=False):

    #     ax.spy(mat)
    from scipy.sparse import issparse
    if issparse(mat):
        im = ax.imshow(mat.todense(), extent=(10, 20, 10, 20), origin='lower', 
                       zorder=1, interpolation='none')

    else:
        im = ax.imshow(mat, extent=(10, 20, 10, 20),
               origin='lower', zorder=1, interpolation='none')


def my_heatmap(mat, fig, ax, title='', vmin=1e-2,vmax=1.0, inverse=True, interactive=False, gene_loc=None,colormap='hot'):

    seq_len = mat.shape[0]

    if interactive is True:
        print ("Interactive, large memory consumer!")
#         plugins.clear(ax)
#         plugins.connect(fig, plugins.Reset(), plugins.BoxZoom(), plugins.Zoom())
        ax.plot(np.arange(seq_len+1, -1), np.arange(seq_len+1, -1),  c='black')
    else:
        ax.plot(np.arange(-1, seq_len+1), np.arange(-1, seq_len+1), c='black')
        if gene_loc is not None:
            print (gene_loc)
            assert len(gene_loc) == 2
            # gene_loc = gene_loc[0]+0.5 + gene_loc[1]+0.5
            ax.plot(np.arange(gene_loc[0]-1, gene_loc[1]), np.arange(gene_loc[0]-1, gene_loc[1]), c='green', linewidth=2)
            ax.plot(np.ones(gene_loc[1]-gene_loc[0]+1)*(gene_loc[1]-1), np.arange(gene_loc[0]-1, gene_loc[1]), c='green', linewidth=1)
            ax.plot(np.arange(gene_loc[0]-1, gene_loc[1]), np.ones(gene_loc[1]-gene_loc[0]+1)*(gene_loc[0]-1),  c='green', linewidth=1)

    if inverse:
        cmap = plt.get_cmap(colormap+'_r')
        if colormap == 'hot':
            cmap = truncate_colormap(cmap, 0.3, 1.1)  # Discards super white range of hit map
    else:
        cmap = plt.get_cmap(colormap)
        if colormap == 'hot':
            cmap = truncate_colormap(cmap, 0.35, 1.0)  # Discards super white range of hit map

#     plugins.connect(fig, plugins.MousePosition(fontsize=14))

    heatmap = ax.matshow(mat, cmap=cmap, vmin=vmin, vmax=vmax, ) #, interpolation='nearest')
#     y, x = np.mgrid[:mat.shape[0], :mat.shape[1]]
#     x,y = x.ravel(),y.ravel()

#     scatter = ax.scatter(x, y, c=mat, s=40, marker='s', edgecolor='none')
#     fig.plugins = [plugins.PointLabelTooltip(scatter, None)]

    ax.set_xticklabels([])
    ax.set_yticklabels([])
    if colormap == 'hot':
        barticks = [vmin] + [r/10.0 for r in range(1, 11)]
        barlabels = [str(vmin)] + [str(r/10.0) for r in range(1, 11)]
    
        cbar = fig.colorbar(heatmap, extend='min', ticks=barticks ,fraction=0.043, pad=0.04)
    else:
        barticks =  [r/10.0 for r in range(-10, 11)]
        barlabels = [str(r/10.0) for r in range(-10, 11)]
        cbar = fig.colorbar(heatmap,  ticks=barticks ,fraction=0.043, pad=0.04)

    # cbar.ax.set_yticklabels(barlabels)
    if inverse:
        cbar.cmap.set_under('white')
    else:
        cbar.cmap.set_under('black')

    #     plt.colorbar(heatmap)

    ax.set_title(title)
    ticks = np.arange(0, mat.shape[0], 10)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    ax.set_xticks(ticks-0.5, minor=True)
    ax.set_yticks(ticks-0.5, minor=True)
#     ax.grid(True, which='minor',color='gray',linewidth=0.0001 )
    ax.grid(False, which='major')  # ,color='gray',linewidth=0.001 )
    #     ax.gca().patch.set_facecolor('0.8')
    ax.tick_params(length=0,
    axis='both',          # changes apply to the x-axis
    which='major',      # both major and minor ticks are affected
    bottom='off',      # ticks along the bottom edge are off
    top='off',         # ticks along the top edge are off
    labelbottom='off')

    ax.set_xlim((-0.5, seq_len-0.5))
#     ax.set_ylim((-0.5,seq_len-0.5))
    ax.set_ylim((seq_len-0.5, -0.5))


def plot_heat_maps_fig(fig, subplot_num, mfe_probs, bp_probs_whole, what='all', inverse=False,
                   interactive=False, gene_loc=None, title_suffix='',vmin=1e-2,vmax=1.0, colormap='hot'):

    if what == 'basepairs' or what == 'all':
        my_heatmap(bp_probs_whole, fig, fig.add_subplot(subplot_num + 1), title_suffix
                   , inverse=inverse, interactive=interactive, gene_loc=gene_loc,vmin=vmin,vmax=vmax,colormap=colormap)

    if what == 'mfe-probs' or what == 'all':
        my_heatmap(mfe_probs, fig, fig.add_subplot(subplot_num + 3), 'mfe:'+title_suffix,
                   inverse=inverse, interactive=interactive, gene_loc=gene_loc,vmin=vmin,vmax=vmax,colormap=colormap)
    if what == 'all':
        my_heatmap(bp_probs_whole*mfe_probs, fig, fig.add_subplot(subplot_num + 2), 'bp*struct:'+title_suffix,
                   inverse=inverse, interactive=interactive, gene_loc=gene_loc,vmin=vmin,vmax=vmax,colormap=colormap)
        my_heatmap(np.sqrt(bp_probs_whole*mfe_probs), fig, fig.add_subplot(subplot_num + 4), 'sqrt(bp*struct):'+title_suffix,
                   inverse=inverse, interactive=interactive, gene_loc=gene_loc,vmin=vmin,vmax=vmax,colormap=colormap)

    #     fig.savefig(filename+'.png', dpi=800)

def plot_heat_maps(mfe_probs, bp_probs_whole, filename='', what='all', inverse=False, 
                   interactive=False, gene_loc=None, title_suffix='',vmin=1e-2,vmax=1.0, out_dir='./',
                   upper_triangle_txt='',lower_triangle_txt='',colormap='hot'):
    if what == 'all':
        fig = plt.figure(figsize=(20, 5))
        subplot_num = 140
    else:
        fig = plt.figure(figsize=(7, 7))
        subplot_num = 110

    if what == 'basepairs' or what == 'all':
        my_heatmap(bp_probs_whole, fig, fig.add_subplot(subplot_num + 1), title_suffix
                   , inverse=inverse, interactive=interactive, gene_loc=gene_loc,vmin=vmin,vmax=vmax,colormap=colormap)

    if what == 'mfe-probs' or what == 'all':
        my_heatmap(mfe_probs, fig, fig.add_subplot(subplot_num + 3), 'mfe-probs:'+title_suffix,
                   inverse=inverse, interactive=interactive, gene_loc=gene_loc,vmin=vmin,vmax=vmax,colormap=colormap)
    if what == 'all':
        my_heatmap(bp_probs_whole*mfe_probs, fig, fig.add_subplot(subplot_num + 2), 'bp*struct:'+title_suffix,
                   inverse=inverse, interactive=interactive, gene_loc=gene_loc,vmin=vmin,vmax=vmax,colormap=colormap)
        my_heatmap(np.sqrt(bp_probs_whole*mfe_probs), fig, fig.add_subplot(subplot_num + 4), 'sqrt(bp*struct):'+title_suffix,
                   inverse=inverse, interactive=interactive, gene_loc=gene_loc,vmin=vmin,vmax=vmax,colormap=colormap)


    fig.text(x=0.72,y=0.3,s=upper_triangle_txt, alpha=0.5)
    fig.text(x=0.65,y=0.18,s=lower_triangle_txt, alpha=0.5)
    
    #     fig.savefig(filename+'.png', dpi=800)
    # if inverse:
        # filename += '_inverse'
    #fig.savefig(filename+'.pdf', dpi=300)

    
    fig.savefig(os.path.join(out_dir, filename+'-dotplot.svg'), dpi=300, format="svg")
    fig.savefig(os.path.join(out_dir, filename+'-dotplot.png'), dpi=600, format="png")
    return fig


def my_heatmaps(rna_seq, context_all, context_len, insert_pos=None, filename='heatmap', what='all', 
                inverse=True, interactive=False, motif_len=None):

    if context_all == None:
        print ("Full sequence and context given")
        whole_seq_context = rna_seq
        gene_loc = [insert_pos, insert_pos+motif_len]
    else:
        context_selection = context_all[0:context_len]
        if insert_pos is None:
            insert_pos = len(context_selection)/2
        whole_seq_context = context_selection[:insert_pos] + rna_seq + context_selection[insert_pos:]
        print (whole_seq_context)
        print (len(rna_seq), len(context_selection), insert_pos)
        gene_loc = [insert_pos, insert_pos+len(rna_seq)]
        
    if what == 'basepairs':
        plot_heat_maps(None, getBPPM(whole_seq_context), filename, what,
                   inverse=inverse, interactive=interactive, gene_loc=gene_loc)
    else:
        plot_heat_maps(get_mfe_probs(whole_seq_context), getBPPM(whole_seq_context), filename, what,
                   inverse=inverse, interactive=interactive, gene_loc=gene_loc)

##################################################################################
# TODO: Merge to parse_dp_ps*( methods it should be easy!
def parse_dp_ps(ps_file):
    '''Extracts base pair probabliies from vienna ps file
    returns: Dictinary of form dict[i:j]=p(i,j) '''

    # Extract sequence from ps file
    myseq = ""
    read_seq = False
    with open(ps_file) as in_ps:
        for line in in_ps:
            if "/sequence" in line:
                read_seq = True
            elif read_seq and ") } def" in line:
                read_seq = False
            elif read_seq:
                myseq += line.rstrip().rstrip("\\")
    #     print ps_file.rstrip("_dp.ps") , myseq

    ureg = re.compile(r'^(\d+)\s+(\d+)\s+(\d+\.\d+)\s+[ul]box\s*')
    bp_prob_dict = dict()
    bp_prob_mat = np.zeros((len(myseq), len(myseq)))

    with open(ps_file) as in_ps:
        for line in in_ps:
            if "ubox" in line or "lbox" in line:
                um = ureg.match(line)
                if um:
                    i, j, sqrp = um.groups()
                    sqrp = float(sqrp)
                    i, j = int(i), int(j)
                    #                     print i, j, sqrp
                    # TODO: Keys and dict not used remov it
                    if "lbox" in line:
                        # keys are pair of indexes as smaller:larger
#                         key = ":".join(sorted([i, j], reverse=False))
                        assert sqrp == 0.95
                        bpprob = 1
                        bp_prob_mat[j-1, i-1] = bpprob
                    else:
                        # keys are pair of indexes as smaller:larger
#                         key = ":".join(sorted([i, j], reverse=True))
                        bpprob = sqrp * sqrp
                        bp_prob_mat[i-1, j-1] = bpprob
#                     assert (key not in bp_prob_dict)
#                     bp_prob_dict[key] = bpprob



    return bp_prob_mat



# Plot ps dotplot files with option to support large files with sparse storage


def parse_dp_ps_sparse(ps_file, sparse=False, bp_range=None, skip_pos=None):
    '''Extracts base pair probabliies from vienna ps file
    returns: Numpy 2d array of form arr[i,j]=p(i,j) '''

    # Extract sequence from ps file
    myseq = ""
    read_seq = False
    with open(ps_file) as in_ps:
        for line in in_ps:
            if "/sequence" in line:
                read_seq = True
            elif read_seq and ") } def" in line:
                read_seq = False
            elif read_seq:
                myseq += line.rstrip().rstrip("\\")
    #     print ps_file.rstrip("_dp.ps") , myseq

    import re

    # e.g. 52 56 0.020043762 ubox
    ureg = re.compile(r'^(\d+)\s+(\d+)\s+(\d+\.\d+)\s+[ul]box\s*')

    # e.g. 1 70 0.9500000 lbox
    #     lreg = re.compile(r'^(\d+)\s+(\d+)\s+(\d+\.\d+)\s+lbox\s*')

    bp_prob_dict = dict()
    mfe_struct_dict = dict()

    if sparse:
        bp_prob_mat = csr_matrix((len(myseq), len(myseq)))

    else:
        bp_prob_mat = np.zeros((len(myseq), len(myseq)))

    with open(ps_file) as in_ps:
        for line in in_ps:
            if "ubox" in line or "lbox" in line:
                um = ureg.match(line)
                if um:
                    i, j, sqrp = um.groups()
                    # keys are pair of indexes as smaller:larger

                    i, j = int(i), int(j)
                    if skip_pos is not None and i >= skip_pos:
                        i += 1
                    if skip_pos is not None and j >= skip_pos:
                        j += 1

                    key = ":".join([str(i), str(j)])

                    if bp_range is not None:  # Skip base pair if one side not in range
                        assert len(bp_range) == 2
                        assert bp_range[0] < bp_range[1]
                        if (i < bp_range[0]-1) or (j > bp_range[1]+1):  # -+1 because not sure for zero based or one based
                            # print i, j, bp_range[0], bp_range[1]
                            continue

                    #                     print i, j, sqrp

                    if "ubox" in line:  # upper triangle of probs
                        assert (key not in bp_prob_dict)
                        bpprob = float(sqrp)*float(sqrp)
                        bp_prob_dict[key] = bpprob

                        bp_prob_mat[i-1, j-1] = bpprob
                    else:  # lower part mfe struct
                        i, j = int(i), int(j)
                        assert (key not in mfe_struct_dict)
                        mfe_struct_dict[key] = 1

    return bp_prob_mat, mfe_struct_dict


def bpp_dict_to_np_array(d, seq):
    np_arr = np.zeros((len(seq), len(seq)))
    for pair in d:
        i, j = pair.split(':')
        np_arr[i, j] = d[pair]
    return np_arr


def plot_dp_ps(dp, sparse=False, gene_loc=None, inverse=False, cut_gene=False, title='', infig=None
               , subplot_num=None):
    from os.path import basename, dirname
#     from pankoff_lib import parse_dp_ps
    
    if cut_gene is True:
        parser_range = gene_loc
    else:
        parser_range = None

    np_arr, mfe_dict = parse_dp_ps_sparse(dp, sparse, parser_range)
    if infig is None:
        fig = plt.figure(figsize=(7, 7))
        subplot_num = 110
    else:
        fig = infig
    if cut_gene is True:
        assert gene_loc is not None
        assert len(gene_loc) == 2
        np_arr_cut = np_arr[gene_loc[0]-2:gene_loc[1]+1, gene_loc[0]-2:gene_loc[1]+1]
        plot_heat_maps_fig(fig, subplot_num, None, np_arr_cut,  what='basepairs', gene_loc=None,
                           inverse=inverse, interactive=False, title_suffix=title)

    else:
        plot_heat_maps_fig(fig, subplot_num, None, np_arr, what='basepairs', gene_loc=gene_loc,
                           inverse=inverse, interactive=False, title_suffix=title)

    if infig is None:
        filename = dirname(dp)+"/"+basename(dp)
        if inverse:
            filename += '_inverse'
#         fig.savefig(filename+'.pdf', dpi=300)
        fig.savefig(filename+'.svg', dpi=300, format="svg")


def plot_differential_dp_ps(dp1, dp2, skips_pos= [None, None], sparse=False, gene_loc=None, inverse=False, cut_gene=False, title='', infig=None
               , subplot_num=None,mut_insert_range=None,mut_delete_range=None,vmin=-1.0,vmax=1.0, ):

    from os.path import basename, dirname
    assert len(skips_pos) == 2
    if cut_gene is True:
        parser_range = gene_loc
    else:
        parser_range = None

#     from pankoff_lib import parse_dp_ps
    np_arr1, mfe_dict1 = parse_dp_ps_sparse(dp1, sparse, parser_range, skips_pos[0])
    np_arr2, mfe_dict2 = parse_dp_ps_sparse(dp2, sparse, parser_range, skips_pos[1])
    if mut_insert_range:
        np_arr2 = np.delete(np.delete(np_arr2,mut_insert_range,0),
          mut_insert_range,1)
    if mut_delete_range:
        np_arr1 = np.delete(np.delete(np_arr1,mut_delete_range,0),
          mut_delete_range,1)


    np_arr = np_arr1 - np_arr2
    if infig is None:
        fig = plt.figure(figsize=(7, 7))
        subplot_num = 110
    else:
        fig = infig
    if cut_gene is True:
        assert gene_loc is not None
        assert len(gene_loc) == 2
        np_arr_cut = np_arr[gene_loc[0]-2:gene_loc[1]+1, gene_loc[0]-2:gene_loc[1]+1]
        plot_heat_maps_fig(fig, subplot_num, None, np_arr_cut,  what='basepairs', gene_loc=None,
                           inverse=inverse, interactive=False, title_suffix=title,vmin=vmin,vmax=vmax,
                           colormap='seismic')

    else:
        plot_heat_maps_fig(fig, subplot_num, None, np_arr, what='basepairs', gene_loc=gene_loc,
                           inverse=inverse, interactive=False, title_suffix=title,vmin=vmin,vmax=vmax,
                           colormap='seismic')

    if infig is None:
        filename = dirname(dp1)+"/"+basename(dp1)+basename(dp2)+'_diff'
        if inverse:
            filename += '_inverse'
#         fig.savefig(filename+'.pdf', dpi=300)
        fig.savefig(filename+'.svg', dpi=300, format="svg")


### Code snippet for Vienna python pavckage
# RNA.cvar.fold_constrained = 0
# structure = ""
# print RNA.pf_fold("AAGGGAAACCCA", structure)
# # RNA.get_subseq_F(4,11)
# f5_p,c_p,fML_p,fM1_p,indx_p,ptype_p = [""]*6
# RNA.export_fold_arrays(f5_p,c_p,fML_p,fM1_p,indx_p,ptype_p)
# # RNA.alloc_sequence_arrays
# from RNA import fold_compound

# comp = RNA.fold_compound
# RNA.export_fold_arrays
# getBPPM("AACGGGAAAGCCA")
