#!/usr/bin/env python

#assign loci preferences

import argparse, os, re, sys
import pandas as pd
import numpy as np
from Bio import SeqIO

def import_negs(infile, loci):
    negl = []
    zerol = []
    df = pd.read_csv(infile)
    for locus in loci:
        filtered_df = df[df['Locus'] == locus]
        if len(filtered_df) == 1:
            negs = filtered_df.iloc[0]['Number_negatives']
            negl.append(negs)
            zeros = filtered_df.iloc[0]['All_zeros']
            zerol.append(zeros)
        else:
            sys.exit(f"{locus}-import_negs: No match or multiple matches found.")
    return negl, zerol


def import_genuscore(infile, loci):
    l = []
    df = pd.read_csv(infile)
    for locus in loci:
        l.append((df.iloc[:, 14:].dropna() == locus).any().any())
    return l

def import_dnds(infile, loci):
    dn = []
    ds = []
    dnds = []
    df = pd.read_csv(infile)
    loci_set = set(df['loci'])
    for locus in loci:
        if locus not in loci_set:
            dn.append(None)
            ds.append(None)
            dnds.append(None)
        else:
            dn.append(df.loc[df['loci'] == locus, 'dn'].iloc[0])
            ds.append(df.loc[df['loci'] == locus, 'ds'].iloc[0])
            dnds.append(df.loc[df['loci'] == locus, 'dnds'].iloc[0])
    return ds, dn, dnds

def import_phastest(infile, loci):
    l = []
    handle = open(infile)
    phastest_loci = []
    search_term = r'EAFFDHGO_[0-9]{5}'
    for line in handle:
        reSearch = re.search(search_term, line)
        if reSearch:
            phastest_loci.append(reSearch.group())
    for locus in loci:
        if locus in phastest_loci:
            l.append(True)
        else:
            l.append(False)
    return l
    

def import_tss(infile, loci):
    l = []
    df = pd.read_csv(infile)
    exported_T3SS = []
    for value in df['NCBI_M12_Xcc']:
        split_value = value.split(', ')
        for locus in split_value:
            exported_T3SS.append(locus)
    for locus in loci:
        if locus in exported_T3SS:
            l.append(True)
        else:
            l.append(False)
    return l

def import_tmhmm(infile, loci):
    l = []
    handle = open(infile)
    tmhmm_loci = []
    search_term = r'EAFFDHGO_[0-9]{5}'
    for line in handle:
        if "Number of predicted TMRs: 0" in line:
            reSearch = re.search(search_term, line)
            if reSearch:
                tmhmm_loci.append(reSearch.group())
    for locus in loci:
        if locus in tmhmm_loci:
            l.append(False)
        else:
            l.append(True)
    return l

def import_signalp(infile, loci):
    l = []
    handle = open(infile)
    signalp_loci = []
    search_term = r'EAFFDHGO_[0-9]{5}'
    for line in handle:
        reSearch = re.search(search_term, line)
        if reSearch:
            signalp_loci.append(reSearch.group())
    for locus in loci:
        if locus in signalp_loci:
            l.append(True)
        else:
            l.append(False)
    return l

def import_biocycle(infile, loci):
    l = []
    df = pd.read_csv(infile)
    biocyc_excluded_loci = df['NCBI_M12_Xcc'].tolist()
    for locus in loci:
        if locus in biocyc_excluded_loci:
            l.append(True)
        else:
            l.append(False)
    return l

def import_cog(infile, loci):
    l =[]
    df = pd.read_csv(infile)
    cog_list = []
    cog_letters = ('X', 'U', 'W', 'V', 'M', 'N')
    for i, value in enumerate(df['COG_LETTER']):
        if value in cog_letters:
            cog_list.append(df.iloc[i]['QUERY_ID'])
    for locus in loci:
        if locus in cog_list:
            l.append(True)
        else:
            l.append(False)
    return l

def import_repeats(infile, loci):
    l = []
    df = pd.read_csv(infile)
    tandem_repeats_loci = df['tandem_repeat'].tolist()
    for locus in loci:
        if locus in tandem_repeats_loci:
            l.append(True)
        else:
            l.append(False)
    return l

def import_homopolymers(seq_d, loci):
    l = []
    for locus in loci:
        s = str(seq_d[locus].seq)
        test = re.findall("(A{8,}|T{8,}|G{8,}|C{8,})",s)
        if len(test) > 0:
            l.append(True)
        else:
            l.append(False)
    return l

def import_gene_sizes(seq_d, loci):
    l = []
    for locus in loci:
        l.append(len(seq_d[locus].seq))
    return l

def get_args():
    parser = argparse.ArgumentParser(description="Generates a True/False matrix and assigns preference to loci from program outputs for filtering of core loci for of MGT scheme generation")

    # Add arguments and options
    parser.add_argument('-n', '--negzeros', required=True, help='File containing number of species negatives and zeros for each loci')
    parser.add_argument('-g', '--genuscore', required=True, help='Gene presence absence file from genus core genome analysis')
    parser.add_argument('-d', '--dnds', required=True, help='Summary file of average dn, ds and dnds values for each loci')
    parser.add_argument('-p', '--phastest', required=True, help='Phastest output')
    parser.add_argument('-t', '--tss', required=True, help='Blast results from Type III secretion system genes')
    parser.add_argument('-m', '--tmhmm', required=True, help='TMHMM output')
    parser.add_argument('-s', '--signalp', required=True, help='SignalP output')
    parser.add_argument('-b', '--biocycle', required=True, help='Biocycle output')
    parser.add_argument('-c', '--cog', required=True, help='COG gene assignments output')
    parser.add_argument('-f', '--fasta', required=True, help='Loci sequences in fasta format')
    parser.add_argument('-r', '--repeat', required=True, help='Repeat finder and homopolyer detection output')
    parser.add_argument('-o', '--outfile', required=True, help='Output prefix for True/False matrix and preference assignments')
    args = parser.parse_args()
    return args

def create_tf_d(args, loci_d):
    d = {}
    d['loci_ids'] = list(loci_d)
    d['negs'], d['zeros'] = import_negs(args.negzeros, d['loci_ids'])
    d['genuscore'] = import_genuscore(args.genuscore, d['loci_ids'])
    d['ds'], d['dn'], d['dnds'] = import_dnds(args.dnds, d['loci_ids'])
    d['phastest'] = import_phastest(args.phastest, d['loci_ids'])
    d['tss'] = import_tss(args.tss, d['loci_ids'])
    d['tmhmm'] = import_tmhmm(args.tmhmm, d['loci_ids'])
    d['signalp'] = import_signalp(args.signalp, d['loci_ids'])
    d['biocycle'] = import_biocycle(args.biocycle, d['loci_ids'])
    d['cog'] = import_cog(args.cog, d['loci_ids'])
    d['gene_size'] = import_gene_sizes(loci_d, d['loci_ids'])
    d['repeat'] = import_repeats(args.repeat, d['loci_ids'])
    d['homopolymer'] = import_homopolymers(loci_d, d['loci_ids'])
    for k, v in d.items():
        print(f'{k}: {len(v)}')
    return d

def generate_combine_boolean_masks(df, d):
    #create boolean masks for determining preferences
    
    mneg0 = df['negs'] == 0
    mneg3 = df['negs'] <= 3
    mneg5 = df['negs'] <= 5
    mzero0 = df['zeros'] == 0
    mzero1 = df['zeros'] <= 1
    mzero2 = df['zeros'] <= 2
    mgcore = df['genuscore'] == True
    mds50 = (df['ds'] >= np.nanpercentile(np.array(d['ds'], dtype=np.float32), 25)) & (df['ds'] <= np.nanpercentile(np.array(d['ds'], dtype=np.float32), 75))
    mds90 = (df['ds'] >= np.nanpercentile(np.array(d['ds'], dtype=np.float32), 5)) & (df['ds'] <= np.nanpercentile(np.array(d['ds'], dtype=np.float32), 95))
    mdnds50 = (df['dnds'] >= np.nanpercentile(np.array(d['dnds'], dtype=np.float32), 25)) & (df['dnds'] <= np.nanpercentile(np.array(d['dnds'], dtype=np.float32), 75))
    mdnds90 = (df['dnds'] >= np.nanpercentile(np.array(d['dnds'], dtype=np.float32), 5)) & (df['dnds'] <= np.nanpercentile(np.array(d['dnds'], dtype=np.float32), 95))
    mphastest = df['phastest'] == False
    mtss = df['tss'] == False
    mtmhmm = df['tmhmm'] == False
    msignalp = df['signalp'] == False
    mbiocycle = df['biocycle'] == False
    mcog = df['cog'] == False
    mgenesize90 = (df['gene_size'] >= np.nanpercentile(np.array(d['gene_size'], dtype=np.float32), 5)) & (df['gene_size'] <= np.nanpercentile(np.array(d['gene_size'], dtype=np.float32), 95))
    mhomopol = df['homopolymer'] == False
    mrepeats = df['repeat'] == False
    print(f"ds 25th Percentile = {np.nanpercentile(np.array(d['ds'], dtype=np.float32), 25)}")
    print(f"ds 75th Percentile = {np.nanpercentile(np.array(d['ds'], dtype=np.float32), 75)}")
    print(f"ds 5th Percentile = {np.nanpercentile(np.array(d['ds'], dtype=np.float32), 5)}")
    print(f"ds 95th Percentile = {np.nanpercentile(np.array(d['ds'], dtype=np.float32), 95)}")
    print(f"dnds 25th Percentile = {np.nanpercentile(np.array(d['dnds'], dtype=np.float32), 25)}")
    print(f"dnds 75th Percentile = {np.nanpercentile(np.array(d['dnds'], dtype=np.float32), 75)}")
    print(f"dnds 5th Percentile = {np.nanpercentile(np.array(d['dnds'], dtype=np.float32), 5)}")
    print(f"dnds 95th Percentile = {np.nanpercentile(np.array(d['dnds'], dtype=np.float32), 95)}")
    print(f"gene_size 5th Percentile = {np.nanpercentile(np.array(d['gene_size'], dtype=np.float32), 5)}")
    print(f"gene_size 95th Percentile = {np.nanpercentile(np.array(d['gene_size'], dtype=np.float32), 95)}")
    #create combined masks for each preference
    masks = {}
    
    masks[1] = mneg0 & mzero0 & mgcore & mds50 & mdnds50 & mphastest & mtss & mtmhmm & msignalp & mbiocycle & mcog & mgenesize90 & mds90 & mdnds90 & mhomopol & mrepeats
    masks[2] = mneg3 & mzero1 & mgcore & mds50 & mdnds50 & mphastest & mtss & mtmhmm & msignalp & mbiocycle & mcog & mgenesize90 & mds90 & mdnds90 & mhomopol & mrepeats
    masks[3] = mneg3 & mzero1 & mds50 & mdnds50 & mphastest & mtss & mtmhmm & msignalp & mbiocycle & mcog & mgenesize90 & mds90 & mdnds90 & mhomopol & mrepeats
    masks[4] = mneg3 & mzero1 & mphastest & mtss & mtmhmm & msignalp & mbiocycle & mcog & mgenesize90 & mds90 & mdnds90 & mhomopol & mrepeats
    masks[5] = mneg3 & mzero1 & mtmhmm & msignalp & mbiocycle & mcog & mgenesize90 & mds90 & mdnds90 & mhomopol & mrepeats
    masks[6] = mneg3 & mzero1 & mbiocycle & mcog & mgenesize90 & mds90 & mdnds90 & mhomopol & mrepeats
    masks[7] = mneg5 & mzero2 & mgenesize90 & mds90 & mdnds90 & mhomopol & mrepeats
    masks[8] = mneg5 & mzero2 & mds90 & mhomopol & mrepeats
    masks[9] = mneg5 & mzero2 & mhomopol & mrepeats
    
    return masks

def assign_loci_to_pref(df, masks):
    prefd_loci = []
    prev_prefs = []
    pref_d = {}
    for pref in range(1, 10):
        pl = df.loc[masks[pref], 'loci_ids'].tolist()
        pref_d[pref] = list(set(pl) - set(prev_prefs))
        prev_prefs += pl
    pref_d[10] = list(set(df['loci_ids']) - set(prev_prefs))
    outpref = {}
    outpref['Genes'] = []
    outpref['Preference'] = []
    for pref in pref_d.keys():
        for Gene in pref_d[pref]:
            outpref['Genes'].append(Gene)
            outpref['Preference'].append(str(pref))
    pref_output_df = pd.DataFrame(outpref)
    return pref_output_df

def main():
    # Add arguments and options
    args = get_args()
    #import reference core gene fasta
    loci_d = SeqIO.to_dict(SeqIO.parse(args.fasta, 'fasta'))
    #import data into combined dictionary
    tf_d = create_tf_d(args, loci_d)
    
    #convert dictionary into dataframe
    tf_df = pd.DataFrame(tf_d)
    masks = generate_combine_boolean_masks(tf_df,tf_d)
    
    #assign loci to preferences
    pref_output_df = assign_loci_to_pref(tf_df, masks)
    
    #output to files
    tf_matrix_file = f"{args.outfile}.result_summary.tsv"
    pref_file = f"{args.outfile}.mgt_loci_preferences.tsv"
    tf_df.to_csv(tf_matrix_file, sep='\t', index=False)
    pref_output_df.to_csv(pref_file, sep='\t', index=False)

if __name__ == "__main__":
    main()
