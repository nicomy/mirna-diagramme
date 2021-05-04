#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 15:45:09 2020

@author: nicolas
"""
import subprocess
from Bio import SeqIO
# from Bio.Seq import Seq
import Bio.Seq as bs
import os 
import numpy as np
import plotnine as pt
import pandas as pd
# import adjustText
# from adjustText import *

import mygene
mg = mygene.MyGeneInfo()

import matplotlib.pyplot as plt

import random
import statistics as st


# tmp_utr_f= "/mnt/HD/Nico/data_PHD/script/output/utr_tmp.fa"    
# tmp_mir_f="/mnt/HD/Nico/data_PHD/script/output/mir_tmp.fa"    
# tmp_out_f= "/mnt/HD/Nico/data_PHD/script/output/out_tmp"   

tmp_utr_f= "/usr/tmp/utr_tmp.fa"    
tmp_mir_f="/usr/tmp/mir_tmp.fa"    
tmp_out_f= "/usr/tmp/out_tmp"   

def create_file(utr_seq,mirna):
    tmp_utr =  open(tmp_utr_f,'w')
    tmp_mir = open(tmp_mir_f,'w')
    tmp_utr.write(">UTR-tmp"+"\n"+utr_seq)
    tmp_mir.write(">mir-mtp\n"+mirna)
    tmp_utr.close()
    tmp_mir.close()
    return(tmp_utr_f,tmp_mir_f,tmp_out_f)

def create_file_tmp(utr_name,mi_name,utr_seq,mirna):
    tmp_utr =  open(tmp_utr_f,'w')
    tmp_mir = open(tmp_mir_f,'w')
    tmp_utr.write(">"+utr_name+"\n"+utr_seq)
    tmp_mir.write(">"+mi_name+"\n"+mirna)
    tmp_utr.close()
    tmp_mir.close()
    return(tmp_utr_f,tmp_mir_f,tmp_out_f)


def inta(utr_seq,mirna,mi,l_start,l_end,l_score=[]):
    # print ("\tcomputing IntaRNA for "+ mi)
    p =     subprocess.Popen(["IntaRNA", "-q",mirna ,"-t",utr_seq,"--out", "STDOUT",'--outMode','C'], 
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = p.communicate(b"input data that is passed to subprocess' stdin")
    l_tmp = output.decode('utf-8').rstrip().split("\n") 
    if (len(l_tmp)<2):
        l_start.append(np.nan)
        l_end.append(np.nan)
        l_score.append(np.nan)
        return()
    l_res = l_tmp[1].split(';')
    if(l_res[1]<l_res[2] ):
        l_start.append( np.int32(l_res[1]))
        l_end.append( np.int32(l_res[2]))
        l_score.append(l_tmp[-1].split(';')[-1])
    else : 
        l_start.append( np.int32(l_res[2]))
        l_end.append( np.int32(l_res[1]))
        l_score.append(l_tmp[-1].split(';')[-1])
    
    # return (l_tmp[1].split(';')[1:3])
        

def miranda(utr_f,mir_f,mi,l_start,l_end,l_score=[]):
    # print ("\tcomputing miranda for "+ mi)
    score=[140,130,120,110,100,90,80,70,60,50,40,30,0]
    miranda_ex = '/usr/Program/miranda/bin/miranda'
    for i in score :     
    #     p =     subprocess.Popen(["miranda", mir_f,utr_f,"-sc",str(i)], 
    #            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p =     subprocess.Popen([miranda_ex, mir_f,utr_f,"-sc",str(i)], 
               stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, err = p.communicate(b"input data that is passed to subprocess' stdin")
        # print(output)
        ss = output.decode('utf-8')
        l_ss= ss.splitlines()
        # print(len(l_ss))
        if (len(l_ss)> 45 ) :
            l_tmp = l_ss[38].rstrip().split(' ')
            if(l_tmp[10][2:]<l_tmp[12]  ):
                l_start.append( np.int32(l_tmp[10][2:]))
                l_end.append( np.int32(l_tmp[12]))
            else : 
                l_start.append( np.int32(l_tmp[12]))
                l_end.append( np.int32(l_tmp[10][2:]))
            l_score.append(l_tmp[4])
            return()
            # return(l_tmp[10][2:],l_tmp[12] ) #score is in l_tmp[4]
            
        
def pita(utr_f,mir_f,mi,out_f,l_start,l_end,l_score=[]):
    
    pita_ex= '/usr/Program/64bit_exe_pita_prediction/pita_prediction.pl'
    p =    subprocess.Popen([pita_ex,'-utr',utr_f,'-mir',mir_f,'-prefix',out_f], 
        stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = p.communicate(b"input data that is passed to subprocess' stdin")
    with open(out_f+"_pita_results.tab") as f:
        lines = f.readlines()
        # print(lines)
        if(len(lines)<2):
            l_start.append(np.nan)
            l_end.append(np.nan)
            l_score.append(np.nan)
            return 

        l_tmp = lines[1].split('\t')
        if(l_tmp[2]<l_tmp[3] ):
            l_start.append( np.int32(l_tmp[2]))
            l_end.append( np.int32(l_tmp[3]))
        else : 
            l_start.append( np.int32(l_tmp[3]))
            l_end.append( np.int32(l_tmp[2]))
        l_score.append(l_tmp[-1].rstrip())

def del_file():
    os.remove(tmp_utr_f)
    os.remove(tmp_mir_f)
    os.remove(tmp_out_f+"_pita_results.tab")
    os.remove(tmp_out_f+"_pita_results_targets.tab")
    os.remove("tmp_seqfile1")
    os.remove("tmp_seqfile2")
    
# =============================================================================
#  utilitaries 
# =============================================================================

def flatten_dic(dic):
    
    def flatten_rec(item,l_tmp):
        if (isinstance(item,str)) :
            l_tmp.append(item)
        else: 
            for sub in item: 
                flatten_rec(sub,l_tmp)
                
                
    for k , v in dic.items() : 
        l_tmp = []
        flatten_rec(v,l_tmp)
        dic[k] = l_tmp 
 
def intersection(lst1, lst2): 
    lst3 = [value for value in lst1 if value in lst2] 
    return lst3     
 
def mklist(n):
    for _ in range(n):
        yield []


def file_size(file,isfasta=True):
    if(isfasta):
        p =     subprocess.Popen(["grep", "-c","^",file], 
                             stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, err = p.communicate(b"input data that is passed to subprocess' stdin")
        total = int(output.decode('utf-8').rstrip())
    else :
        
        p =     subprocess.Popen(["grep", "-c","^",file], 
                         stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, err = p.communicate(b"input data that is passed to subprocess' stdin")
        total = int(output.decode('utf-8').rstrip())
    return (total)


def randomize_interation(total,nb_item,randomize=True):
    if(randomize):
        set_item = set()
        while(len(set_item) < nb_item):
             set_item.add(random.randint(1,total))
        tmp = list(set_item)
        tmp.sort()
        return(tmp)     
    else : 
        return(list(range(1,nb_item+1)))

def clean_list(l_tmp):
    l_res = []
    for i in l_tmp :
         if (i):
             if('\t' in i ):
                 j  = i.split('\t')
                 for x in j :
                     if(x):
                         l_res.append(x.rstrip())
             else : 
                 l_res.append(i)
    return(l_res)


def remove_nan(l):
    res=[]
    for i in l : 
        if (not np.isnan(i)):
            res.append(i)
    return res

# create a motif from the alternate letter : ex base_motif="ACBAC" , letter =" A C  " res = "AABCC"
def add_motif(base_motif, letter):
    if (len(base_motif)!= len(letter)):
        raise ("motif len not equal")
    base_motif = list(base_motif)
    for i in range(0,len(base_motif)):
        if(letter[i] != ' '):
            base_motif[i]=letter[i] 
    return ''.join(base_motif)

# =============================================================================
# stat computation 
# =============================================================================

def stat(l):
    l_tmp = remove_nan(l)
    
    if(l_tmp):
        mean = st.mean(l_tmp)
        v_min = min(l_tmp)
        v_max= max(l_tmp)
        variance = st.pstdev(l_tmp)
        return([mean,v_min,v_max,variance])
    return ([np.nan, np.nan,np.nan,np.nan])
# =============================================================================
# Function to add to graphic
# =============================================================================

adjust_text_dict = {
    'expand_points': (2, 2),
    'arrowprops': {
        'arrowstyle': '->',
        'color': 'red'
    }
}

def add_mirna_g(g,df, str_name,str_start,str_end,dis_pos,l_s,l_e,l_score=[]):
    # print(str_name,str_start,str_end,dis_pos,l_s,l_e)
    df[str_start]= pd.Series(l_s)
    df[str_end] = pd.Series(l_e)

    g+= pt.annotate("text", x=0,y=dis_pos,label=str_name)
    g+= pt.geom_errorbarh(df,pt.aes(xmin=str_start,y=(dis_pos),xmax=str_end,color='mi_name'))
    g+= pt.geom_segment(df,pt.aes(x=str_start,y=(dis_pos),yend=0,xend=str_start,color='mi_name'))
    if(l_score):
        # print(l_score)
        # pd.options.display.float_format = '{:.1f}'.format
        score_column_name = 'score'+str_name
        # print(l_score,score_column_name,str_start,dis_pos)
        df[score_column_name] = pd.Series(l_score,dtype=np.float).map('{:.0f}'.format)
        
        g+= pt.geom_text(df, pt.aes(x=str_start,y=dis_pos,label=score_column_name,color='mi_name'),
                          nudge_x=0.1, nudge_y=0.1)#,adjust_text=adjust_text_dict)
        # g+= pt.geom_text(df, pt.aes(x=str_start,y=dis_pos,label=score_column_name,color='mi_name'))

    
# =============================================================================
#  Function to create dictionnary and stuff 
# =============================================================================



def miRNA_dictionnary(miRNA_file):
    d_mir = {}
    tot_size_mir =0 
    mir_fasta= SeqIO.parse(open(miRNA_file),'fasta') 
    for mir  in mir_fasta: 
        tmp = str(mir.seq)
        d_mir[mir.id] = tmp
        tot_size_mir += len(tmp)
    return(d_mir)        

def open_dic(dic_file):
    dic={}
    with open(dic_file) as f:  
        for l in f : 
           l_tmp = l.split("\t")     
           tmp  = l_tmp[1].split(",")
           if (tmp[-1]):
               dic[l_tmp[0]] = tmp[:-1]
               dic[l_tmp[0]].append(tmp[-1].rstrip())
           else : 
               dic[l_tmp[0]] = tmp[:-1]
    return(dic)

def Te2Gid(gid2Te_file):
    d_te2gid={}
    with open(gid2Te_file) as f:  
        for gid in f : 
           l_tmp = gid.split("\t")     
           tmp  = l_tmp[1].split(",")
           if (tmp[-1]):
               l_TE = tmp[:-1]
               l_TE.append(tmp[-1].rstrip())
           else : 
               l_TE = tmp[:-1]
           if(l_TE):
           # d_gid2Te[l_tmp[0]] = l_TE
               for te in l_TE : 
                   d_te2gid[te] = l_tmp[0] 
           else :
               d_te2gid[l_tmp[1].rstrip()] = l_tmp[0] 
    return(d_te2gid)
       
    
def gid2mir(miRNA2gid,d_mir) :
    d_gid2mir = {}
    with open(miRNA2gid) as f:  
        for mir in f :
            l_tmp = mir.split(" ")
            # d_mi2gid[l_tmp[0]] = l_tmp[1].rstrip()
            mi_name = l_tmp[0]
            if(mi_name in d_mir):
                gid =  l_tmp[1].rstrip()
                if (gid in d_gid2mir):
                    d_gid2mir[gid].append(l_tmp[0])
                else :
                    d_gid2mir[gid]=[l_tmp[0]]
            else : 
                continue 
    return(d_gid2mir)

def write_d(output_file, d, seq = False):
    with open(output_file,'w') as out :
        for k,v in d.items() :
            if (seq):
                out.write(">"+k+"\n"+v+"\n")
            else : 
                out.write(str(k))
                out.write("\t")
                if(isinstance(v, list)):
                    for i in v[:-1] : 
                        out.write(str(i))
                        out.write(",")
                    out.write(str(v[-1]))
                else:
                    out.write(str(v))
                out.write("\n")


def get_seed(dic_mir):
    dic_seed = {}
    for mi, seq in dic_mir.items():
        dic_seed[mi]= bs.Seq(seq[1:8]).reverse_complement()
    return (dic_seed)

def get_smile_output(smile_output1):
    d_smile1={}
    
    smile1 = open(smile_output1,'r') 
    line = smile1.readline()
    while("================" not in line):
        line = smile1.readline()
    line = smile1.readline()
    
    while ('Nb models' not in line) :
        l_tmp = line.split(' ')
        l_tmp2 = l_tmp[-1].split('\t')
        
        d_smile1[l_tmp[0]] = (l_tmp2[0],l_tmp2[1].rstrip()) 
        line = smile1.readline()
    smile1.close()
    return (d_smile1)

def get_smile_output_val(smile_output2):
    d_smile_val1occ = {}
    d_smile_val = {}
    smile1 = open(smile_output2,'r') 
    line = smile1.readline()
    while("================" not in line):
        line = smile1.readline()
    line = smile1.readline()
    
    while (line != "\n") :
        l_tmp = clean_list(line.split(' '))
        d_smile_val1occ[l_tmp[0]]=l_tmp[1:]
        line = smile1.readline()
    
    while("================" not in line):
        line = smile1.readline()
    line = smile1.readline()
    while ("User time" not in line ) :
        l_tmp = clean_list(line.split(' '))
        d_smile_val[l_tmp[0]]=l_tmp[1:]
        line = smile1.readline()

    return(d_smile_val1occ,d_smile_val)

# =============================================================================
# functions to request ncbi database
# =============================================================================


def getEnsemble_ID(l_gid,d_TR,d_gE,d_TE,ensembl_missing)  : 
    count= 0 
    total = len(l_gid)    
    refseq_missing = []
    missing_data = [] 
    
    for g_id in l_gid:
        count +=1 
        print(count,"/",total)
        # g_id = int(gid.strip())
        t =mg.getgene(g_id,'refseq,ensembl')
        if(t):
            if('ensembl' in t):
                if(isinstance(t['ensembl'],dict)):
                    d_gE[g_id] = t['ensembl']['gene']
                    if(isinstance(t['ensembl']['transcript'],list)):
                        d_TE[g_id] = t['ensembl']['transcript']
                    else : 
                        d_TE.setdefault(g_id,[]).append(t['ensembl']['transcript'])
                else : 
                    # g_id are uniq so wre do not need to check if they are already in the dic
                    d_gE[g_id]=[]
                    d_TE[g_id]=[]
                    for e in t['ensembl'] :
                        d_gE[g_id].append(e['gene'])
                        d_TE[g_id].append(e['transcript'])
            else :
                ensembl_missing.append(g_id)
                
            if('refseq' in t ):
                if('rna' in t['refseq']):
                    d_TR[g_id] = t['refseq']['rna']
                else: 
                    refseq_missing.append(g_id)
            else : 
                refseq_missing.append(g_id)
        else : 
            missing_data.append(g_id)
            
    flatten_dic(d_TE)
    print("missing data -> Ensembl missing: "+
          str(len(ensembl_missing)) +" \t Refseq missing: "+ str(len(refseq_missing)) )

# =============================================================================
# function for box _plot 
# =============================================================================

def box_plot(list_to_plot,title):
    plt.boxplot(list_to_plot,notch ='True')
    # ax.set_yticklabels(cnt_name[0:5])
    plt.title(title +" with outlier")
    plt.show()
    
    plt.boxplot(list_to_plot,notch ='True',showfliers=False)
    # ax.set_yticklabels("")
    plt.title(title+" without outlier")
    plt.show()