#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
***********************************************************
* Author  : Zhou Wei                                      *
* Date    : 2020/07/27 09:58:45                           *
* E-mail  : welljoea@gmail.com                            *
* Version : --                                            *
* You are using the program scripted by Zhou Wei.         *
* If you find some bugs, please                           *
* Please let me know and acknowledge in your publication. *
* Thank you!                                              *
* Best wishes!                                            *
***********************************************************
'''
import argparse
import os
import re
import traceback
import numpy as  np
import pandas as pd
import glob
from joblib import Parallel, delayed
import seaborn as sns
import matplotlib.pyplot as plt

class Getinfo():
    def __init__(self,args):
        self.args = args
        self.IDs = [ os.path.basename(i) for i in glob.glob(self.args.input + '/*') if os.path.exists(i + '/SS2/STAR') ]

class Mapping(Getinfo):
    def starstate(self, Sid):
        starfile = '%s/%s/SS2/STAR/%s__D1__L1.Log.final.out' %(self.args.input, Sid, Sid)
        try:
            finfo = open(starfile, 'r').readlines()[5:]
            finfo = [ [i.strip(), '-'] if len(i.split('\t'))==1  
                        else re.split( ' \|\t', i.strip())
                        for i in finfo]
            finfo = pd.DataFrame(finfo, columns=['mappingstate', Sid])
            finfo.set_index('mappingstate',inplace=True)
            return finfo
        except IOError:
            print('%s is not exit.')
    def starmerge(self):
        Starpd = Parallel( n_jobs=-1)( delayed( self.starstate )( _n ) for _n in self.IDs )
        Starpd = pd.concat(Starpd, axis=1).T
        Starpd.to_csv(self.args.outdir + '/mappingstate.xls', sep='\t',header=True, index=True, index_label='sampleID')


class RSEMce(Getinfo):
    def unitpm(self, Sid, redup='Max'):
        rsemgene = '%s/%s/SS2/RSEM/%s__D1.rsemce.genes.results' %(self.args.input, Sid, Sid)
        rsemgene = pd.read_csv(rsemgene,header=0,sep='\t')
        rsemgene.gene_id = rsemgene.gene_id.str.split('_').str[-1]
        rsemgene = rsemgene[(rsemgene.TPM >0)][['gene_id','TPM']].copy()

        if redup=='Sum':
            rsemgene=rsemgene.groupby(by='gene_id', sort=False).sum()
        elif redup=='Max':
            rsemgene=rsemgene.groupby(by='gene_id', sort=False).max()

        rsemgene.rename(columns={'TPM':Sid+'_TPM'}, inplace=True)
        return rsemgene
    def mergetpm(self):
        rsemtpm = Parallel( n_jobs=-1)( delayed( self.unitpm )( _n ) for _n in self.IDs )
        rsemtpm = pd.concat(rsemtpm, axis=1,sort=False)
        rsemtpm.to_csv(self.args.outdir + '/rsem.tpm.xls', sep='\t',header=True, index=True, index_label='Gene')

    def tpmplot(self):
        rsemtpm  = pd.read_csv(self.args.outdir + '/rsem.tpm.xls', header=0, sep='\t', index_col=0 )
        countsdf = pd.melt(rsemtpm.reset_index(),
                             id_vars='Gene',
                             value_vars=rsemtpm.columns.to_list(),
                             var_name='sampleid',
                             value_name='TPM')
        countsdf = countsdf[countsdf.TPM >0]
        plt.figure(figsize=(15,7))
        vs = sns.violinplot(x="sampleid", y="TPM", data=countsdf)
        #vs = sns.swarmplot(x="sampleid", y="TPM", data=countsdf,edgecolor="none", s=2, color="black")
        vs.set_yscale('log')
        plt.xticks( rotation='270')
        plt.savefig( self.args.outdir + '/rsem.tpm.counts.pdf',  bbox_inches='tight' )
        plt.close()

        rsemtpm1 = rsemtpm.copy()
        rsemtpm1[rsemtpm1<1]=np.nan
        genecount = pd.concat([rsemtpm.count(axis=0), rsemtpm1.count(axis=0) ],axis=1 )
        genecount.reset_index(inplace=True)
        genecount.columns=['sampleID','genecount','genecountlg1']
        genecount.to_csv(self.args.outdir + '/rsem.tpm.genecounts.xls', sep='\t',header=True, index=False)

        if os.path.exists(self.args.outdir + '/mappingstate.xls'):
            Starpd=pd.read_csv(self.args.outdir + '/mappingstate.xls', header=0, sep='\t', index_col=0 )
            tpmcou=genecount.copy()
            tpmcou.sampleID=tpmcou.sampleID.str.strip("_TPM")
            tpmcou.set_index('sampleID',inplace=True)
            QC_Gene=pd.concat( (Starpd, tpmcou), axis=1 )
            QC_Gene.to_csv(self.args.outdir + '/mapping.QC.Gcounts.xls', sep='\t',header=True, index=True, index_label='sampleID')

        genecount = genecount.melt(id_vars=['sampleID'], 
                             value_vars=['genecount','genecountlg1'],
                             var_name='count_type',
                             value_name='counts')
        genecount.sort_values(['sampleID'], inplace=True)
        plt.figure(figsize=(12,8))
        bt = sns.barplot(x='sampleID',
                            y='counts',
                            hue='count_type',
                            hue_order=['genecount','genecountlg1'], 
                            palette=sns.color_palette("Set2"), data=genecount)
        plt.xticks( rotation='270')
        leg = plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), numpoints=1)
        plt.savefig( self.args.outdir + '/rsem.tpm.genecounts.pdf',bbox_extra_artists=(leg,), bbox_inches='tight' )
        plt.close()


def Args():
    parser = argparse.ArgumentParser(
                formatter_class=argparse.RawDescriptionHelpFormatter,
                prefix_chars='-+',
                conflict_handler='resolve',
                description="\nmapping QC state :\n",
                epilog='''Example:''')

    parser.add_argument("-i", "--input",     type=str,  required=True, help="the input file including Id,snv.vcf.gz and bam columns")
    parser.add_argument("-o", "--outdir",    type=str,  default=os.getcwd(), help="output file dir, default=current dir")
    parser.add_argument("-p", "--pool",      type=int,  default=23,  help="the CPU numbers that can be used")
    parser.add_argument("-a", "--analysis",  type=str,  default="annovar,Lichee,sciclone,pyclone", help="the analysis type, default=Lichee,choices=annovar,Lichee,sciclone,pyclone")
    parser.add_argument("-d", "--do",        action="store_true",default=False, help='do the workshell')
    args  = parser.parse_args()
    return args
def Commands():
    args = Args()
    os.makedirs( os.path.dirname(args.outdir) , exist_ok=True)
    print("The argument you have set as follows:".center(59, '*'))
    for i,k in enumerate(vars(args),start=1):
        print('**%s|%-13s: %s'%(str(i).zfill(2), k, str(getattr(args, k))) )
    print(59 * '*')

    try:
        Mapping(args).starmerge()
        RSEMce(args).mergetpm()
        RSEMce(args).tpmplot()
        print('Success!!!')
    except Exception:
        traceback.print_exc()
    finally:
        print('You can check your progress in log file.')
Commands()
