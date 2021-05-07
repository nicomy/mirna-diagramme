#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 11:07:07 2021

@author: nhomberg
"""

import numpy as np
import pandas as pd
import RnaFun as RF
import itertools

import plotly.graph_objects as go
from ipywidgets import widgets
import Bio.Seq as bs
    
def read_clash_csv2df(f_clash_data_csv) :
    df_clash = pd.read_csv(f_clash_data_csv,skiprows=30,sep="\t")
    df_clash['nb_targets'] = df_clash.groupby(['microRNA_name'])["mRNA_name"].transform('count')
    return(df_clash)
    
    
# =============================================================================
#  reading smile : 

    # Para : 
    #     smile_f_z, 
    #     df_clash, 
    #     min_len,
    #     max_len,
    #     min_target
    
    # return: 
    #     df : smile in df
# =============================================================================
    
def parce_read_smile(file_smile_z, df_clash, min_len,max_len,min_target):
    # d_smile1 = RF.get_smile_output(smile_output1)

    
    # smile_name_val_1occ = [ '%right','right'	, '%shfl.','shfl.','Sigma','Chi2',	'Z-score']
    # smile_name_val =['right','shfl.','Sigma','Chi2','Z-score']
    
    d_smile_val1occ, d_smile_val = RF.get_smile_output_val(file_smile_z)
    
    # columns_name = ['mirna','Seed','#seed in UTR 1occ','#seed in UTR','#mi interactions',
                    # "Chi2 1occ","Z_score 1occ","chi2","Z_score"]
    
    columns_name_df = ['right','shfl.','Sigma','Chi2','Z_score','mirna','#interactions']
    
    
# =============================================================================
# crete set from sub-string of the mirna
# =============================================================================

    d_submot = {}
    
    
    # df_targets = df_clash[df_clash['nb_targets']>=min_target][['microRNA_name','nb_targets']].drop_duplicates().set_index(["microRNA_name"])
    df_targets = df_clash[['microRNA_name','nb_targets']].drop_duplicates().set_index(["microRNA_name"])
    # start_time = timeit.default_timer()
    # df_targets.at['MIMAT0017994_MirBase_miR-3615_microRNA','nb_targets']
    
    for row in df_clash[["microRNA_name","miRNA_seq"]].drop_duplicates().itertuples():
        _,mi_name,seq = row
        # print(mi)
        seq = str(bs.Seq(seq).reverse_complement())
        for len_i in range(min_len,max_len+1):
            # print(len_i)
            len_seq = len(seq)
            for i in range(0,len_seq-len_i+1):
                # print("\t"+seq[i:len_i+i])
                mot = seq[i:len_i+i]
                # name_short= mi_name.split('_')[2]
                # d_submot[ mot ]  = d_submot.get(mot, []) + [mi_name]
                d_submot.setdefault(mot, []).append(mi_name)
    
    
    # count = 0 
    # for k,v in d_submot.items():
    #     if len(v) > 1 : 
    #         count +=1
    # print(count,len(d_submot),sep='\t' )
    # compare smile with submotif inside miRna   
    nb_not_found =  0
    s_mirna_found = set()
    # tmp_name = 'MIMAT0004518_MirBase_miR-16-2*_microRNA'
    for motif, l_v in d_smile_val.items() :
        if (motif in d_submot):
            l_mir = d_submot[motif]
            d_smile_val[motif].append(l_mir)
            s_mirna_found.update(l_mir)
            l_inter = []
            for mir in l_mir : 
                l_inter.append(df_targets.at[mir,'nb_targets'])
            d_smile_val[motif].append(l_inter)
        else : 
            d_smile_val[motif].append(np.nan)
            nb_not_found+=1
    
    # =============================================================================
    # create smile data frame linked with mirna and display 
    # =============================================================================
    df = pd.DataFrame.from_dict(d_smile_val, orient ='index',columns=columns_name_df)
    df['Z_score'] = pd.to_numeric(df['Z_score'])
    df.sort_values(['Z_score'],ascending=1)
    pos = np.arange(1,len(df)+1) 
    df['position']=pos
    # df.to_csv(csv_output)
    dt = df.dropna()
    
    s_mrna_all_found = set()
    for l_mi in dt["mirna"]:
        s_mrna_all_found.update(l_mi)
    
    # number_of_useless_motif_all = df['mirna'].isnull().sum()
    return(df)


# df_smile_short_seq[((df_smile_short_seq['Z_score' ]> thres_zscore_short) 
#                           |  (df_smile_short_seq['Z_score' ] < -thres_zscore_short))
#                           & df_smile_short_seq['#interactions'].notnull()]

def select_Zscore(df_smile,thres_h_zscore=None, thres_l_zscore=None,clean_nmatch=True):
    equation = pd.Series([True]*len(df_smile),index=df_smile.index)
    if (clean_nmatch):
        equation = df_smile['#interactions'].notnull()
    if(thres_h_zscore ):
        if(thres_l_zscore) :
            equation = equation & ((df_smile['Z_score' ]> thres_h_zscore ) | (df_smile['Z_score' ] < -thres_l_zscore))
        else :
            equation = equation & (df_smile['Z_score' ]> thres_h_zscore )
    elif(thres_l_zscore) :
            equation = equation & (df_smile['Z_score' ] < -thres_l_zscore)

    return( df_smile[equation]) 


def mirna_set(df):
    return (set(itertools.chain(*df['mirna'].dropna())))
    
                      
def create_db_for_bar_plotly(df1,df2,name1,name2,Zscore1,Zscore2):
    
    df1_all = select_Zscore(df1,Zscore1,Zscore1,False)
    df2_all = select_Zscore(df2,Zscore2,Zscore2,False)
    df_inter= pd.merge(df1_all, df2_all,left_index=True, right_index=True, suffixes=('', '_drop')) 
    df_inter.drop([col for col in df_inter.columns if 'drop' in col], axis=1, inplace=True)
    
    df1_match = select_Zscore(df1_all)
    df2_match = select_Zscore(df2_all)
    # df_inter_match = pd.merge(df1_match, df2_match,left_index=True, right_index=True, suffixes=('', '_drop')) 
    df_inter_match= select_Zscore(df_inter)
    # df_inter_match.drop([col for col in df_inter_match.columns if 'drop' in col], axis=1, inplace=True)
    
    
    l_name=[name1,"intersection",name2]
    
    count_all = [len(df1_all),len(df_inter),len(df2_all)]
    # count_match = [len(df1_all['#interactions'].notnull()),len(df_inter['#interactions'].notnull()),len(df2_all['#interactions'].notnull())]
    count_match = [len(df1_match),len(df_inter_match),len(df2_match)]
    nb_mirna= [len(mirna_set(df1_all)),len(mirna_set(df1_all)),len(mirna_set(df2_all))]
    
    dic = {"all":count_all,'matching':count_match,"name":l_name,"nb_mirna":nb_mirna}
    return( pd.DataFrame.from_dict(dic).set_index('name'))
    

def response(change,g,df_smile_short_seq,df_smile_fullseq,use_zscore_short,use_zscore_long,thres_zscore_short,thres_zscore_long):
    if (use_zscore_long.value):
        z_score_long = thres_zscore_long.value
    else : 
        z_score_long = None
    if (use_zscore_short.value):
        z_score_short = thres_zscore_short.value 
    else : 
        z_score_short = None
    
    df_plot=create_db_for_bar_plotly(df_smile_short_seq,df_smile_fullseq,
                                        "short_seq","long_seq",z_score_short,z_score_long)

    with g.batch_update():
        g.data[0].y = df_plot['all']
        g.data[0].customdata =   np.transpose([df_plot['nb_mirna']])   
        g.data[0].texttemplate="%{y}  mirnas : %{customdata[0]}"
        
        g.data[1].y = df_plot['matching']
        g.data[1].customdata = np.transpose([(df_plot['matching']/ df_plot['all'])*100 ]) 
        g.data[1].texttemplate="%{y} (%{customdata[0]:.2f}%)"
        

def plotly_bar_with(f_seq_short,f_seq_longue,clash_csv,motif_min,motif_max):
    
    min_target= 10     
    thres_zscore_short_value =10 
    thres_zsocre_long_value = 20 
    
    # =============================================================================
    # import de donnée 
    # =============================================================================
    
    df_clash = read_clash_csv2df(clash_csv)
    
    nb_mirna_total = len(set(df_clash.microRNA_name))
    
    
    df_smile_short_seq = parce_read_smile(f_seq_short,df_clash,motif_min,motif_max,min_target)
    df_smile_fullseq = parce_read_smile(f_seq_longue,df_clash,motif_min,motif_max,min_target)
    
    
    #non nécessaire
    # dfb_short = select_Zscore(df_smile_short_seq, thres_zscore_short_value,thres_zscore_short_value)
    # dfb_long = select_Zscore(df_smile_fullseq,thres_zsocre_long_value,thres_zsocre_long_value)
    
    # intersected_df = pd.merge(dfb_short, dfb_long,left_index=True, right_index=True)    
    
    pd.options.plotting.backend = "plotly"


    df_plot=create_db_for_bar_plotly(df_smile_short_seq,df_smile_fullseq,
                                        "short_seq","long_seq",thres_zscore_short_value,thres_zsocre_long_value)
    
    
    #print(df_plot)
    
    #df = pd.DataFrame(dict(a=[1,3,2], b=[3,2,1]))
    #df2 = pd.DataFrame(dict(a=[10,30,20], b=[30,20,10]))
    
     
    fig1 = go.Bar(x= df_plot.index, y = df_plot['all'],name='all motif',
                  customdata =  np.transpose([df_plot['nb_mirna']]),        
                  texttemplate="%{y}  mirnas : %{customdata[0]}  ",
                  hovertemplate="<br>".join([
                "total: %{y}","nombre de mirna trouvé :%{customdata[0]} "
            ]),textposition='auto')
    fig2 = go.Bar(x= df_plot.index, y = df_plot['matching'],name ='matching motif',
                  customdata = np.transpose([(df_plot['matching']/ df_plot['all'])*100 ]),
                  texttemplate="%{y}   (%{customdata[0]:.2f}%)",
                textposition='auto')
    #fig3 = go.Bar(x= df_plot.index, y = df_plot['matching'])
    
    
    
    # slicer et case
    thres_zscore_short = widgets.IntSlider(
        value=10.0,
        min=0.0,
        max=50.0,
        step=1.0,
        description='Z-score limite :',
        continuous_update=False
    )
    use_zscore_short = widgets.Checkbox(
        description='filter z-score for short_seq: ',
        value=True,
    )
    
    container = widgets.HBox(children=[use_zscore_short, thres_zscore_short])
    
    thres_zscore_long = widgets.IntSlider(
        value=20.0,
        min=0.0,
        max=60.0,
        step=1.0,
        description='Z-score limite long :',
        continuous_update=False
    )
    use_zscore_long = widgets.Checkbox(
        description='filter z-score for long_seq: ',
        value=True,
    )
    
    container2 = widgets.HBox(children=[use_zscore_long, thres_zscore_long])
    
    
    
    
    g = go.FigureWidget(
                        data=[fig1,fig2],
                        layout=go.Layout(
                            title=dict(
                                text='Motif founds on short and long séquences, and the number of mirnas found /%{nb_mirna_total} mirnas total'
                            ),
                            barmode='overlay'
                        ))
    
    
    l_v_box = [  container,container2, g]
    l_var = [use_zscore_short,use_zscore_long,thres_zscore_short,use_zscore_long]
    return (l_var,l_v_box)
    # fonction d'update du graphique 

            
            
    # thres_zscore_short.observe(response, names="value")
    # use_zscore_short.observe(response, names="value")
    # thres_zscore_long.observe(response, names="value")
    # use_zscore_long.observe(response, names="value")
    
    # widgets.VBox([  container,container2, g]) 
        
    
    
    
    
                      
