from time import sleep as sl
import random, sys, os
import pandas as pd
import argparse
from Bio import SeqIO

random.seed(145156156)

def too_close(poslists,newpos,limit):
    """

    :param poslists: list of start positions
    :param newpos: start position of gene to test
    :param limit: minimum distance cutoff
    :return: True if within limit, false if further
    """
    newpos = int(newpos)
    for i in poslists:
        if -1*limit <= i-newpos <= limit:
            return True
    return False

def import_starts(genbank_file):
    starts = {}
    recs = SeqIO.parse(genbank_file, 'gb')
    for rec in recs:
        for feat in rec.features:
            if feat.type == 'CDS':
                starts[feat.qualifiers['locus_tag'][0]] = int(feat.location.start)
    return starts

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Assign genes to MGT levels based on preference input')
    parser.add_argument('--input_path', type=str, required=True, help='Path to the preference input file.')
    parser.add_argument('--manually_assigned_list_file', type=str, required=True, help='Path to the list of manually assigned MGT genes.')
    parser.add_argument('--output_folder', type=str, required=True, help='Path to the output folder. Will create if it doesn\'t exist')
    parser.add_argument('--genbank_file', type=str, required=True, help='Path to the MGT reference sequence genbank file.')
    args = parser.parse_args()
    
    filt = pd.read_csv(args.input_path)
    manually_assigned_list = args.manually_assigned_list_file
    outfolder = args.output_folder
    starts = import_starts(args.genbank_file)

os.makedirs(outfolder, exist_ok=True)
filt = filt.merge(pd.DataFrame.from_dict(starts, orient='index', columns=['start']), left_on='loci_ids', right_index=True)

# #scheme target sizes
target_sizes = {'MGT4':129578,'MGT5':259157,'MGT6':1036626}

#scheme lowest allowed loci preference numbers
preflimit = {'MGT4':4,'MGT5':6,'MGT6':8}

#scheme smallest distance allowed between loci
distlimit = {'MGT4':10000,'MGT5':2000,'MGT6':500}

#read in manually assigned genes
assigned_genes = {}
with open(manually_assigned_list) as mal_h:
    assigned_genes['MGT3'] = []
    for line in mal_h:
        assigned_genes['MGT3'].append(line.strip())

outputs = assigned_genes

#list of used genes
donegenes = []
for key, value in assigned_genes.items():
    for i in value:
        donegenes.append(i)

#length of total input genes
input_length = {}
for key, value in assigned_genes.items():
    length = 0
    for el in value:
        #print(el, key, value, filt[filt['loci_ids'] == el]['gene_size'])
        #print(filt)
        gene_length = int(filt[filt['loci_ids'] == el]['gene_size'].values)
        length = length + gene_length
    input_length[key] = length

#list of gene start pos
startposs = []
for key, value in assigned_genes.items():
    for el in value:
        start = int(filt[filt['loci_ids'] == el]['start'].values)
        startposs.append(start)

#dict of genes to preferences
prefassigns = {}
for el in donegenes:
    prefassigns[el] = 1

#dict of distance to list of loci that are less that that distance
toclose_genes = {x:[] for x in range(0,20001,1000)}

#preference number
prefno = 1
#length limit
limno = 20000

#get initial list of loci names from preference 1
pref_loci = list(filt[filt["preference"]==prefno]["loci_ids"])

#make a dataframe of oly preference 1
pref_df = pd.DataFrame(filt[filt["preference"]==prefno])

for i in sorted(target_sizes.keys()):
    print(i)
    if not i in list(outputs):
        outputs[i] = []
        totlen = 0
        locino = 0
    else:
        outputs[i] = assigned_genes[i]
        totlen = input_length[i]
        locino = len(assigned_genes[i])
    # outputs[i] = []
    # totlen = 0
    # locino = 0

    #while current running total of scheme genes is less than the target size
    while totlen < target_sizes[i]:
        #if there are any genes left in the pref_loci list
        if len(pref_loci) > 0:
            #get number of genes
            geneno = len(pref_loci)
            #get a random number from 0 to above number
            a = random.randint(0, geneno-1)

            #get gene name
            ids = pref_loci[a]
            #get gene info from df
            gene = filt.loc[filt["loci_ids"]==ids]
            pos = int(gene["start"])
            length = int(gene["gene_size"])
            # if chro not in startposs:
            #     startposs[chro] = []
            #check if new gene is too close to existing picks given current limit (limno) using start position of new gene
            to_close_filt = too_close(startposs,pos,limno)
            ##if the new gene is not already picked and is not already noted as too close fr current limit
            if ids not in donegenes and ids not in toclose_genes[limno]:
                #if the new gene is not too close
                if to_close_filt == False:
                    #add gene to done list, outputs for current MGT, start position list, preference recording, total length so far for MGT, number of loci, remove gene from current preference list
                    donegenes.append(ids)
                    outputs[i].append(ids)
                    startposs.append(pos)
                    prefassigns[ids] = prefno
                    totlen += length
                    locino+=1
                    pref_loci.remove(ids)
                else:
                    ## if gene is too close add to to_close_genes dictionary for current limit
                    toclose_genes[limno].append(ids)
                    ## remove gene from current preference list
                    pref_loci.remove(ids)
            else:
                #if gene is already in done or toclose remove from current preference list
                pref_loci.remove(ids)
        #if there are no more genes in current preference and the preference number is under the limit for the MGT
        elif prefno < preflimit[i]:
            # print("HELLLOOOO2")
            #increase the preference no by 1
            prefno +=1
            # get loci for this preference
            pref_loci = list(filt[filt["preference"]==prefno]["loci_ids"])
            pref_df = pd.DataFrame(filt[filt["preference"] == prefno])
            # print("preference_no: "+str(prefno))
        # if length limit is  is more than cutoff for min distance
        elif limno > distlimit[i]:
            # reset preference to 1 and reset pref_loci to preference 1
            prefno = 1
            pref_loci = list(filt[filt["preference"]==prefno]["loci_ids"])
            #reduce distance limit to 1000 less
            limno -= 1000
            #if limit is more than 0
            if limno > 0:
                #add lower cuttoff genes to current - because genes in 6000 cutoff chould be included in 7000 cutoff etc
                toclose_genes[limno] += toclose_genes[limno-1000]

            # print(toclose_genes)
            print("distance limit: "+str(limno),len(toclose_genes[limno]))
        else:
            break
    #test for getting loci when not enough sequence is included
    if totlen < target_sizes[i]:
        print("NOT ENOUGH SEQ!!\nNeed: "+str(target_sizes[i]))
        print(i,totlen,locino)

        sl(100000)
    #reset distance limit for each scheme
    limno = 20000


outsummary = open(outfolder + "/all_schemes_loci.txt","w")
out_alleles = os.path.join(outfolder, 'refonly_allelic_profiles')
os.makedirs(out_alleles, exist_ok=True)
for i in outputs:
    outf = open(outfolder +"/"+i+"_gene_accessions.txt","w")
    print(i,len(outputs[i]))
    for gene in outputs[i]:
        outsummary.write("{}\t{}\t{}\n".format(i,gene,prefassigns[gene]))
    outf.write("\n".join(outputs[i]))
    outf.close()
outsummary.close()

#writeout allelic profiles
for i in outputs:
    with open(out_alleles +"/"+i+"_gene_profiles.txt","w") as out:
        out.write('ST' + '\t' + 'dST' + '\t')
        for gene in outputs[i]:
            out.write(gene + '\t')
        out.write('\n')
        out.write('1' + '\t' + '0' + '\t')
        for gene in outputs[i]:
            out.write('1' + '\t')

