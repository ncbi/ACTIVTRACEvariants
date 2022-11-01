# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 14:50:54 2021

@author: connorrp
for ILLUMINA-ONT test set analysis
"""


def typer(array):

    if len(array["ref"]) != len(array["alt"]) or\
        (array["ref"] == "-") or (array["alt"] == "-"):
        return "InDel"
    elif len(array["ref"]) == len(array["alt"]) and\
        (array["ref"] != "-") and (array["alt"] != "-"):
        return "SNP"
    else:
        print(array)
        raise


def combiner():

    import pandas as pd
    import os

    os.chdir("Cleaned_Trace_ONT")

    NCBI = pd.read_csv("NCBI_results7.csv",
                       sep=",",
                       header=0,
                       index_col=False)
    BEI = pd.read_csv("BEI_results8.csv",
                      sep=",",
                      header=0,
                      index_col=False)
    BEI.Acc = BEI.Acc.apply(lambda x: x.split('.')[0])
    Galaxy = pd.read_csv("Galaxy_results3.csv",
                         sep=",",
                         header=0,
                         index_col=False)
    Gilead = pd.read_csv("Gilead_results6.csv",
                         sep=",",
                         header=0,
                         index_col=False)
    LANL = pd.read_csv("LANL_results4.csv",
                       sep=",",
                       header=0,
                       index_col=False)
    LANL.alt = LANL.alt.apply(lambda x: x.split("(")[0])
    LANL = LANL[LANL.ref != LANL.alt]
    Vir = pd.read_csv("VIR_results2.csv",
                      sep=",",
                      header=0,
                      index_col=False)
    Lilly = pd.read_csv("Lilly_results2.csv",
                        sep=",",
                        header=0,
                        index_col=False)
    Map = pd.read_csv("acc_map.csv",
                      sep=",",
                      header=0,
                      index_col=False)

    df = pd.concat([NCBI, BEI, Gilead, LANL,
                    Vir, Lilly, Galaxy],
                   ignore_index=True,
                   axis=0)
    df = df[(df.AF >= 0.15) & (df.Acc.isin(Map.acc))]
    
    df["var"] = df.Acc+"_"+df.ref+df.pos.astype(str)+df.alt
    df["type"] = df.apply(lambda x: typer(x), axis=1)
    platform_dict = Map[['acc',
                         'platform']].\
                         drop_duplicates().set_index('acc').to_dict()
    df["platform"] = df["Acc"].apply(lambda x: platform_dict['platform'][x])
    sample_dict = Map[['acc',
                       'biosample']].\
                       drop_duplicates().set_index('acc').to_dict()
    df["biosample"] = df["Acc"].apply(lambda x: sample_dict['biosample'][x])
    df = df[(df['pos']>500) & (df['pos']<29500)]
    
    software = {'aligner': {'NCBI': 'HISAT2', 'BEI': 'BWA-MEM',
                            'Galaxy':'BWA-MEM', 'Gilead': 'SMALT',
                            'LANL': 'BWA-MEM', 'VIR': 'BWA-MEM',
                            'Lilly': 'BWA-MEM'},
                'caller': {'NCBI': 'GATK', 'BEI': 'FreeBayes',
                           'Galaxy': 'LoFreq', 'Gilead': 'custom',
                           'LANL': 'samtools', 'VIR': 'LoFreq',
                           'Lilly': 'FreeBayes'}}
    df['Aligner'] = df.Group.apply(lambda x: software['aligner'][x])
    df['Caller'] = df.Group.apply(lambda x: software['caller'][x])

    
    print("combined", df.index.size - df.dropna().index.size)
    print("map", Map.index.size - Map.dropna().index.size)
    
    shared_acc = set(NCBI.Acc).intersection(
                    set(BEI.Acc).intersection(
                    set(Galaxy.Acc).intersection(
                    set(Gilead.Acc).intersection(
                    set(Lilly.Acc).intersection(
                    set(LANL.Acc).intersection(
                    set(Vir.Acc)))))))
    shared_acc2 = set(NCBI.Acc).intersection(
                    set(BEI.Acc).intersection(
                    set(Galaxy.Acc).intersection(
                    set(Lilly.Acc).intersection(
                    set(Gilead.Acc).intersection(
                    set(LANL.Acc))))))
    print(df.Acc.unique().size, len(shared_acc), len(shared_acc2))
    return df.dropna(), Map.dropna()#, shared_acc, shared_acc2


def counter(df1, df2, fname):

    import pandas as pd

    df1 = df1[["Group", "var", "type", "platform"]]
    df3 = pd.DataFrame(columns=["Platform", "Groups", "Type", "Count"])

    for p in df2.platform.unique():
        df4 = df1[df1.platform == p]
        for t in df4["type"].unique():
            df5 = df4[df4["type"] == t]
            for v in df5["var"].unique():
                groups = "_".join(set(
                        df5[(df5["var"] == v)].Group
                        ))
                try:
                    idx = df3[(df3.Type == t)
                              & (df3.Groups == groups)
                              & (df3.Platform == p)].index.item()
                    df3.at[idx, "Count"] += 1
                except:
                    df3 = df3.append({"Platform": p,
                                      "Groups": groups,
                                      "Type": t,
                                      "Count": 1},
                                     ignore_index=True)

    print("counter", df3.index.size - df3.dropna().index.size)
    _fname = ''.join(['groupings', fname, '.csv'])
    df3.to_csv(_fname)
    return df3.dropna()


def counter2(df1):

    import pandas as pd

    idx = []
    for bs in df1.biosample.unique():
        for g in df1.Group.unique():
            df = df1[(df1.Group == g) & (df1.biosample == bs)]
            if len(df.platform.unique()) != 2:
                idx = idx + df.index.to_list()
    df1 = df1[~(df1.index.isin(idx))]

    df1 = df1[["Group", "var", "type", "platform"]]
    df3 = pd.DataFrame(columns=["Group", "Type", "Platform", "Count"])

    for g in df1.Group.unique():
        df4 = df1[df1.Group == g]
        for t in df4["type"].unique():
            df5 = df4[df4["type"] == t]
            vars_ = []
            for v in df5["var"].unique():
                var = v.split("_")[1]
                if var not in vars_:
                    df6 = df5[df5["var"].str.contains(var)]
                    if len(df6.platform.unique()) == 2:
                        platform = "Both"
                    elif len(df6.platform.unique()) == 1:
                        platform = df6.platform.unique().item()
                    try:
                        idx = df3[(df3.Type == t)
                                  & (df3.Group == g)
                                  & (df3.Platform == platform)].index.item()
                        df3.at[idx, "Count"] += 1
                    except:
                        df3 = df3.append({"Group": g,
                                          "Type": t,
                                          "Platform": platform,
                                          "Count": 1},
                                         ignore_index=True)
                    vars_.append(var)

    print("counter2", df3.index.size - df3.dropna().index.size)
    return df3.dropna()


def counter3(df1, df2, fname):

    import pandas as pd

    df1 = df1[["Acc", "Groups", "var", "type", "platform",
               'Group', 'Aligner', 'Caller', 'Aligners', 'Callers']]
    df3 = pd.DataFrame(columns=["Platform", "Groups", "Type", "Count"])
    dff3 = pd.DataFrame(columns=["Platform", "Aligners", "Type", "Count"])
    dfff3 = pd.DataFrame(columns=["Platform", "Callers", "Type", "Count"])
    df1["var"] = df1["var"].apply(lambda x: x.split('_')[1])

    for p in df2.platform.unique():
        df4 = df1[df1.platform == p]
        for t in df4["type"].unique():
            df5 = df4[df4["type"] == t]
            for v in df5["var"].unique():
                df6 = df5[df5["var"] == v]
                v_tot = df6.Acc.unique().size
                for g in df6.Groups.unique():
                    df7 = df6[df6.Groups == g]
                    v_sum = df7.Acc.unique().size/v_tot
                    try:
                        idx = df3[(df3.Type == t)
                                & (df3.Groups == g)
                                & (df3.Platform == p)].index.item()
                        df3.at[idx, "Count"] += v_sum
                    except:
                        df3 = df3.append({"Platform": p,
                                        "Groups": g,
                                        "Type": t,
                                        "Count": v_sum},
                                        ignore_index=True)
                for f, d in [('Aligners', dff3), ('Callers', dfff3)]:
                    for s in df6[f].unique():
                        df7 = df6[df6[f] == s]
                        g_tot = df7[df7[f[:-1]].isin(s.split('_'))].\
                                    Group.unique().size
                        v_sum = (df7.Acc.unique().size/g_tot)/v_tot
                        #print(g_tot, v_sum)
                        if f == 'Aligners':
                            try:
                                idx = dff3[(dff3.Type == t)
                                        & (dff3.Aligners == s)
                                        & (dff3.Platform == p)].index.item()
                                dff3.at[idx, "Count"] += v_sum
                                #print(dff3.index.size)
                            except:
                                dff3 = dff3.append({"Platform": p,
                                              "Aligners": s,
                                              "Type": t,
                                              "Count": v_sum},
                                                ignore_index=True)
                                #print(dff3.index.size)
                        elif f == 'Callers':
                            try:
                                idx = dfff3[(dfff3.Type == t)
                                        & (dfff3.Callers == s)
                                        & (dfff3.Platform == p)].index.item()
                                dfff3.at[idx, "Count"] += v_sum
                                #print(d.index.size)
                            except:
                                dfff3 = dfff3.append({"Platform": p,
                                              "Callers": s,
                                              "Type": t,
                                              "Count": v_sum},
                                                ignore_index=True)
                                #print(dfff3.index.size)

    print("counter", df3.index.size - df3.dropna().index.size)
    _fname = ''.join(['groupings_by_var', fname, '.csv'])
    df3.to_csv(_fname)
    _fname = ''.join(['alinger_groupings_by_var', fname, '.csv'])
    dff3.to_csv(_fname)
    _fname = ''.join(['caller_groupings_by_var', fname, '.csv'])
    dfff3.to_csv(_fname)
    return df3.dropna(), dff3.dropna(), dfff3.dropna()


def grouper(df1, df2):

    from statistics import mean
    df1.reset_index(inplace=True, drop=True)
    for p in df2.platform.unique():
        df_ = df1[df1.platform == p]
        for t in df_['type'].unique():
            df__ = df_[df_['type']==t]
            for v in df__['var'].unique():
                groups = '_'.join(
                        set(df__[df__['var']==v]['Group'])
                        )
                avg_AF = mean(df__[(df__["var"] == v)].AF)
                avg_DP = mean(df__[(df__["var"] == v)].DP)
                idx = df__[df__['var']==v].index
                for i in idx:
                    df1.at[i, "Groups"] = groups
                    df1.at[i, "avg_AF"] = avg_AF
                    df1.at[i, "avg_DP"] = avg_DP
    print("grouper", df1.index.size - df1.dropna().index.size)

    vars_ = []
    for t in df1['type'].unique():
        df4 = df1[df1['type'] == t]
        for v in df4["var"].unique():
            v = v.split('_')[1]
            if v not in vars_:
                df5 = df4[df4["var"].str.contains(v)]
                for g in df5.Group.unique():
                    platforms = "_".join(set(
                            df5[(df5.Group == g)].platform
                            ))
                    idx = df1[(df1["type"] == t)
                              & (df1["var"].str.contains(v))
                              & (df1.Group == g)].index
                    for i in idx:
                        df1.at[i, "Platforms"] = platforms
                vars_.append(v)
                
    for f in ['Aligner', 'Caller']:
            vars_ = []
            for t in df1['type'].unique():
                df4 = df1[(df1['type'] == t)]
                for v in df4["var"].unique():
                    v = v.split('_')[1]
                    if v not in vars_:
                        df5 = df4[df4["var"].str.contains(v)]
                        software = list(set(df5[f]))
                        software.sort()
                        software = "_".join(software)
                        idx = df1[(df1["type"] == t)
                                    & (df1["var"].str.contains(v))].index
                        for i in idx:
                            df1.at[i, f+'s'] = software
                        vars_.append(v)

    print("grouper", df1.index.size - df1.dropna().index.size)
    return df1.dropna()


def boxplotter(df, fname):

    from matplotlib import pyplot as plt
    import seaborn as sns

    print(df.columns)
    df.reset_index(inplace=True, drop=True)
    print(df.columns)
    s = df.Groups.str.len().sort_values().index

    data_ = df.reindex(s)

    plt.figure(figsize=(20, 10))
    chart = sns.boxplot(x="Groups",
                        y="avg_AF",
                        data=data_[data_["type"] == "SNP"],
                        hue="platform")
    plt.xticks(rotation=45, horizontalalignment='right')
    plt.tight_layout()
    _fname1 = ''.join(['snp_hist', fname, '.png'])
    plt.savefig(_fname1)
    plt.clf()

    plt.figure(figsize=(20, 10))
    chart = sns.boxplot(x="Groups",
                        y="avg_AF",
                        data=data_[data_["type"] == "InDel"],
                        hue="platform")
    plt.xticks(rotation=45, horizontalalignment='right')
    plt.tight_layout()
    _fname2 = ''.join(['indel_hist', fname, '.png'])
    plt.savefig(_fname2)
    plt.clf()
    del plt, chart

def upsetplotter(df, fname):

    from matplotlib import pyplot as plt
    from upsetplot import from_memberships
    from upsetplot import plot

    groupings = [x.split("_") for x in
                 df[(df["Platform"] == "ILLUMINA")
                 & (df["Type"] == "SNP")].Groups]
    data = from_memberships(groupings,
                            data=[x for x in
                                  df[(df["Platform"] == "ILLUMINA")
                                  & (df["Type"] == "SNP")]["Count"]])
    plt.figure(figsize=(20, 10))
    fig = plot(data)
    plt.tight_layout()
    _fname1 = ''.join(['ILLUMINDA_SNP_upset', fname, '.png'])
    plt.savefig(_fname1)
    plt.clf()

    groupings = [x.split("_") for x in
                 df[(df["Platform"] == "ILLUMINA")
                 & (df["Type"] == "InDel")].Groups]
    data = from_memberships(groupings,
                            data=[x for x in
                                  df[(df["Platform"] == "ILLUMINA")
                                  & (df["Type"] == "InDel")]["Count"]])
    plt.figure(figsize=(20, 10))
    fig = plot(data)
    plt.tight_layout()
    _fname2 = ''.join(['ILLUMINA_InDel_upset', fname, '.png'])
    plt.savefig(_fname2)
    plt.clf()

    groupings = [x.split("_") for x in
                 df[(df["Platform"] == "OXFORD_NANOPORE")
                 & (df["Type"] == "SNP")].Groups]
    data = from_memberships(groupings,
                            data=[x for x in
                                  df[(df["Platform"] == "OXFORD_NANOPORE")
                                  & (df["Type"] == "SNP")]["Count"]])
    plt.figure(figsize=(20, 10))
    fig = plot(data)
    plt.tight_layout()
    _fname3 = ''.join(['ONT_SNP_upset', fname, '.png'])
    plt.savefig(_fname3)
    plt.clf()

    groupings = [x.split("_") for x in
                 df[(df["Platform"] == "OXFORD_NANOPORE")
                 & (df["Type"] == "InDel")].Groups]
    data = from_memberships(groupings,
                            data=[x for x in
                                  df[(df["Platform"] == "OXFORD_NANOPORE")
                                  & (df["Type"] == "InDel")]["Count"]])
    plt.figure(figsize=(20, 10))
    fig = plot(data)
    plt.tight_layout()
    _fname4 = ''.join(['ONT_InDel_upset', fname, '.png'])
    plt.savefig(_fname4)
    plt.clf()
    del plt, fig


def barplotter(df, fname):

    from matplotlib import pyplot as plt
    import seaborn as sns
    import matplotlib.patches as mpatches

    for t in df['Type'].unique():
        for g in df['Group'].unique():
            total = df[(df["Type"] == t) & (df["Group"] == g)]["Count"].sum()
            for i in df[(df["Type"] == t) & (df["Group"] == g)].index:
                df.at[i, 'Count'] = df.at[i, 'Count']/total
                #print(df.at[i, 'Count'])
    plt.figure(figsize=(20, 10))
    chart1 = sns.barplot(x="Group",
                        y="Count",
                        data=df[df["Type"] == "SNP"].\
                                groupby('Group')['Count'].sum().reset_index(),
                        color='darkblue')
    chart2 = sns.barplot(x="Group",
                        y="Count",
                        data=df[(df["Type"] == "SNP") &
                                (df.Platform != "OXFORD_NANOPORE")].\
                                groupby('Group')['Count'].sum().reset_index(),
                        color='blue')
    chart3 = sns.barplot(x="Group",
                        y="Count",
                        data=df[(df["Type"] == "SNP") &
                                (df.Platform == "Both")].\
                                groupby('Group')['Count'].sum().reset_index(),
                        color='lightblue')
    top_bar = mpatches.Patch(color='darkblue', label='ONT')
    middle_bar = mpatches.Patch(color='blue', label='Illumina')
    bottom_bar = mpatches.Patch(color='lightblue', label='Both')
    plt.legend(handles=[top_bar, middle_bar, bottom_bar])
    plt.xticks(rotation=45, horizontalalignment='right')
    plt.tight_layout()
    _fname1 = ''.join(['SNP_platform_compare', fname, '.png'])
    plt.savefig(_fname1)
    plt.clf()

    plt.figure(figsize=(20, 10))
    chart1 = sns.barplot(x="Group",
                        y="Count",
                        data=df[df["Type"] == "InDel"].\
                                groupby('Group')['Count'].sum().reset_index(),
                        color='darkblue')
    chart2 = sns.barplot(x="Group",
                        y="Count",
                        data=df[(df["Type"] == "InDel") &
                                (df.Platform != "OXFORD_NANOPORE")].\
                                groupby('Group')['Count'].sum().reset_index(),
                        color='blue')
    chart3 = sns.barplot(x="Group",
                        y="Count",
                        data=df[(df["Type"] == "InDel") &
                                (df.Platform == "Both")].\
                                groupby('Group')['Count'].sum().reset_index(),
                        color='lightblue')
    top_bar = mpatches.Patch(color='darkblue', label='ONT')
    middle_bar = mpatches.Patch(color='blue', label='Illumina')
    bottom_bar = mpatches.Patch(color='lightblue', label='Both')
    plt.legend(handles=[top_bar, middle_bar, bottom_bar])
    plt.xticks(rotation=45, horizontalalignment='right')
    plt.tight_layout()
    _fname2 = ''.join(['InDel_platform_compare', fname, '.png'])
    plt.savefig(_fname2)
    plt.clf()
    del plt, chart1, chart2, chart3, top_bar, middle_bar, bottom_bar


def scatterplotter(df, fname):

    from matplotlib import pyplot as plt
    import seaborn as sns

    plt.figure(figsize=(20, 10))
    df2 = df[~df.Group.isin(['Vir', 'Lilly'])]
    df2 = df2[df2['type'] == 'SNP']
    df2.Groups = df2.Groups.apply(lambda x: len(x.split('_')))
    df2 = df2[['var', 'Groups', 'Platforms',
               'avg_AF', 'avg_DP']].drop_duplicates()
    chart = sns.scatterplot(data=df2,
                            x="avg_DP",
                            y="avg_AF",
                            hue="Groups",
                            style="Platforms",
                            palette="viridis")
    chart.set(xscale='log')
    plt.tight_layout()
    _fname1 = ''.join(['SNP_scatter', fname, '.png'])
    plt.savefig(_fname1)
    plt.clf()

    plt.figure(figsize=(20, 10))
    df2 = df[~df.Group.isin(['Vir', 'Lilly'])]
    df2 = df2[df2['type'] == 'SNP']
    df2.Groups = df2.Groups.apply(lambda x: len(x.split('_')))
    df2 = df2[df2.Groups >= 3]
    df2 = df2[['var', 'Groups', 'Platforms',
               'avg_AF', 'avg_DP']].drop_duplicates()
    chart = sns.scatterplot(data=df2,
                            x="avg_DP",
                            y="avg_AF",
                            hue="Groups",
                            style="Platforms",
                            palette="viridis")
    chart.set(xscale='log')
    plt.tight_layout()
    _fname2 = ''.join(['SNP_scatter_majority', fname, '.png'])
    plt.savefig(_fname2)
    plt.clf()

    plt.figure(figsize=(20, 10))
    df2 = df[~df.Group.isin(['Vir', 'Lilly'])]
    df2 = df2[df2['type'] == 'SNP']
    df2.Groups = df2.Groups.apply(lambda x: len(x.split('_')))
    df2 = df2[['var', 'platform', 'Groups', 'Platforms',
               'avg_AF', 'avg_DP']].drop_duplicates()
    chart = sns.relplot(data=df2,
                        x="avg_DP",
                        y="avg_AF",
                        col="platform",
                        hue="Groups",
                        style="Platforms",
                        palette="viridis",
                        kind='scatter')
    chart.set(xscale='log')
    try:
        chart._legend.set_bbox_to_anchor(1)
    except:
        pass
    plt.tight_layout()
    _fname3 = ''.join(['scatter_byplat', fname, '.png'])
    plt.savefig(_fname3)
    plt.clf()

    plt.figure(figsize=(20, 10))
    df2 = df[~df.Group.isin(['Vir', 'Lilly'])]
    df2 = df2[df2['type'] == 'SNP']
    df2.Groups = df2.Groups.apply(lambda x: len(x.split('_')))
    df2 = df2[df2.Groups >= 3]
    df2 = df2[['var', 'platform', 'Groups',
               'Platforms', 'avg_AF', 'avg_DP']].drop_duplicates()
    chart = sns.relplot(data=df2,
                        x="avg_DP",
                        y="avg_AF",
                        col="platform",
                        hue="Groups",
                        style="Platforms",
                        palette="viridis",
                        kind='scatter')
    chart.set(xscale='log')
    try:
        chart._legend.set_bbox_to_anchor(1)
    except:
        pass
    plt.tight_layout()
    _fname4 = ''.join(['scatter_byplat_majority', fname, '.png'])
    plt.savefig(_fname4)
    plt.clf()
    del plt, chart
    
    
def clustermapper(df5, fname):
    
    import numpy as np
    import pandas as pd
    from matplotlib import pyplot as plt
    import seaborn as sns

    df7 = pd.DataFrame(columns=
                       ["platform", "pos", "type", "BEI",
                        "LANL", "NCBI", "VIR", "Lilly", 
                        "Gilead", "Galaxy", "var"])
    for p in df5.platform.unique():
        df8 = df5[df5.platform == p]
        for t in df8.type.unique():
            df9 = df8[df8.type == t]
            for v in df9['var'].unique():
                df10 = df9[df9['var'] == v]
                v = v.split('_')[1]
                bei, lanl, ncbi, vir,\
                gilead, galaxy, lilly =\
                    0, 0, 0, 0, 0, 0, 0
                groups = set([x for y in df10.Groups.unique()
                             for x in y.split('_')])
                pos = df10.pos.unique().item()
                if "BEI" in groups:
                    bei = 1
                if "LANL" in groups:
                    lanl = 1
                if "NCBI" in groups:
                    ncbi = 1
                if "VIR" in groups:
                    vir = 1
                if "Lilly" in groups:
                    lilly = 1
                if "Galaxy" in groups:
                    galaxy = 1
                if "Gilead" in groups:
                    gilead = 1
                try:
                    idx = df7[(df7.type == t)
                              & (df7['var'] == v)
                              & (df7.platform == p)
                              & (df7.pos == pos)].index.item()
                    df7.at[idx, "BEI"] += bei
                    df7.at[idx, "LANL"] += lanl
                    df7.at[idx, "NCBI"] += ncbi
                    df7.at[idx, "VIR"] += vir
                    df7.at[idx, "Lilly"] += lilly
                    df7.at[idx, "Galaxy"] += galaxy
                    df7.at[idx, "Gilead"] += gilead
                except:
                    df7 = df7.append({"platform": p,
                                      "pos": pos,
                                      "var": v,
                                      "type": t,
                                      "BEI": bei,
                                      "LANL": lanl,
                                      "NCBI": ncbi,
                                      "VIR": vir,
                                      "Lilly": lilly,
                                      "Galaxy": galaxy,
                                      "Gilead": gilead},
                                     ignore_index=True)

    plt.figure(figsize=(120, 80))
    data_ = df7[(df7.platform == "ILLUMINA")
                & (df7.type == 'SNP')][["pos", "var", "BEI", "LANL",
                                        "NCBI", "VIR", "Lilly",
                                        "Gilead", "Galaxy"]]
    data_.pop('var')
    data_.pos = data_.pos.astype(int)
    data_ = data_.set_index("pos").sort_index()
    data_ = data_ + np.random.normal(0,
                                     1e-10,
                                     [data_.index.size,data_.columns.size])
    g = sns.clustermap(data_.dropna().astype(float),
                       cmap='vlag', metric="euclidean", z_score=0,
                       row_cluster=False, linewidths=0)
    plt.tight_layout()
    _fname1 = ''.join(['Illumina_SNP_heatmap', fname, '.png'])
    plt.savefig(_fname1, bbox_inches='tight', dpi=300)
    plt.clf()

    plt.figure(figsize=(120, 80))
    data_ = df7[(df7.platform == "ILLUMINA")
                & (df7.type == 'InDel')][["pos", "var", "BEI", "LANL",
                                          "NCBI", "VIR", "Lilly",
                                          "Gilead", "Galaxy"]]
    data_.pop('var')
    data_.pos = data_.pos.astype(int)
    data_ = data_.set_index("pos").sort_index()
    data_ = data_ + np.random.normal(0,
                                     1e-10,
                                     [data_.index.size,data_.columns.size])
    g = sns.clustermap(data_.dropna().astype(float),
                       cmap='vlag', metric="euclidean", z_score=0,
                       row_cluster=False, linewidths=0)
    plt.tight_layout()
    _fname2 = ''.join(['Illumina_InDel_heatmap', fname, '.png'])
    plt.savefig(_fname2, bbox_inches='tight', dpi=300)
    plt.clf()

    plt.figure(figsize=(120, 80))
    data_ = df7[(df7.platform == "OXFORD_NANOPORE")
                & (df7.type == 'SNP')][["pos", "var", "BEI", "LANL",
                                        "NCBI", "Gilead", "Galaxy"]]
    data_.pop('var')
    data_.pos = data_.pos.astype(int)
    data_ = data_.set_index("pos").sort_index()
    data_ = data_ + np.random.normal(0,
                                     1e-10,
                                     [data_.index.size,data_.columns.size])
    g = sns.clustermap(data_.dropna().astype(float),
                       cmap='vlag', metric="euclidean", z_score=0,
                       row_cluster=False, linewidths=0)
    plt.tight_layout()
    _fname3 = ''.join(['ONT_SNP_heatmap', fname, '.png'])
    plt.savefig(_fname3, bbox_inches='tight', dpi=300)
    plt.clf()

    plt.figure(figsize=(120, 180))
    data_ = df7[(df7.platform == "OXFORD_NANOPORE")
                & (df7.type == 'InDel')][["pos", "var", "BEI", "LANL",
                                          "Gilead", "Galaxy"]]
    data_.pop('var')
    data_.pos = data_.pos.astype(int)
    data_ = data_.set_index("pos").sort_index()
    data_ = data_ + np.random.normal(0,
                                     1e-10,
                                     [data_.index.size,data_.columns.size])
    g = sns.clustermap(data_.dropna().astype(float),
                       cmap='vlag', metric="euclidean", z_score=0,
                       row_cluster=False, linewidths=0)
    plt.tight_layout()
    _fname4 = ''.join(['ONT_InDel_heatmap', fname, '.png'])
    plt.savefig(_fname4, bbox_inches='tight', dpi=300)
    plt.clf()
    del plt, g
    
    
def scorer(df, fname):
    
    import pandas as pd
    
    dfout = pd.DataFrame(columns = ['group', 'platform',
                                    'type', 'total', 'fraction'])
    for p in df.platform.unique():
        for t in df['type'].unique():
            df2 = df[(df.platform == p) & (df['type'] == t)]
            tote = df2['var'].unique().size
            n_groups = df2.Groups.unique().size
            majority = 0
            for v in df2['var'].unique():
                if len(df2[df2['var'] == v].\
                       Groups.unique().item().split('_')) >= n_groups/2:
                    majority += 1
            i = dfout.index.size
            dfout.at[i, 'group'] = 'overall'
            dfout.at[i, 'platform'] = p
            dfout.at[i, 'type'] = t
            dfout.at[i, 'total'] = majority
            dfout.at[i, 'fraction'] = majority/tote
            dfout.at[i, 'len'] = None
            for g in df2['Group'].unique():
                count = df2[df2['Group'] == g]['var'].unique().size
                percent = count/tote
                i = dfout.index.size
                dfout.at[i, 'group'] = g
                dfout.at[i, 'platform'] = p
                dfout.at[i, 'type'] = t
                dfout.at[i, 'total'] = count
                dfout.at[i, 'fraction'] = percent
                dfout.at[i, 'len'] = None
            if t == 'InDel':
                df2['idlen'] = 0
                for i in df2.index:
                    diff = abs(len(df2.at[i, 'ref']) - len(df2.at[i, 'alt']))
                    if df2.at[i, 'alt'] == "-":
                        df2.at[i, 'idlen'] = diff + 1
                    else:
                        df2.at[i, 'idlen'] = diff
                for l in df2.idlen.unique():
                    df3 = df2[df2.idlen==l]
                    tote = df3['var'].unique().size
                    n_groups = df3.Groups.unique().size
                    majority = 0
                    for v in df3['var'].unique():
                        if len(df3[df3['var'] == v].\
                               Groups.unique().item().split('_')) \
                               >= n_groups/2:
                            majority += 1
                    i = dfout.index.size
                    dfout.at[i, 'group'] = 'overall'
                    dfout.at[i, 'platform'] = p
                    dfout.at[i, 'type'] = t
                    dfout.at[i, 'total'] = majority
                    dfout.at[i, 'fraction'] = majority/tote
                    dfout.at[i, 'len'] = l
                    for g in df3['Group'].unique():
                        count = df3[df3['Group'] == g]['var'].unique().size
                        percent = count/tote
                        i = dfout.index.size
                        dfout.at[i, 'group'] = g
                        dfout.at[i, 'platform'] = p
                        dfout.at[i, 'type'] = t
                        dfout.at[i, 'total'] = count
                        dfout.at[i, 'fraction'] = percent
                        dfout.at[i, 'len'] = l
                
    dfout.to_csv(fname+'_score.csv')


def reprocess(inpt, df2, fname, AF=[0], DP=[0], ALTDP=[0],
              HQ=[False], HP=[None], SE=[None], BC=[None]):
    
    import pandas as pd
    
    Map = df2
    
    for af in AF:
        for dp in DP:
            for altdp in ALTDP:
                for hq in HQ:
                    for hp in HP:
                        for se in SE:
                            for bc in BC:
                                fname0 = fname.split('_')[0]
                                fname2 = '_'.join(fname.split('_')[1:])
                                fname1 = 'AF'+str(int(100*af))\
                                +'DP'+str(dp)\
                                +'ALTDP'+str(altdp)
                                if hq:
                                    fname1 = fname1+"HQ"
                                if hp in [True, False]:
                                    if hp:
                                        fname1 = fname1+"HP"
                                    else:
                                        fname1 = fname1+"NHP"
                                if se in [True, False]:
                                    if se:
                                        fname1 = fname1+"SE"
                                    else:
                                        fname1 = fname1+"PE"
                                if bc in ['pos', 'bs']:
                                    if bc == 'pos':
                                        fname1 = fname1+"pos"
                                    if bc == 'bs':
                                        fname1 = fname1+"bs"
                                fname_ = fname0+fname1+'_'+fname2
                                print(fname, '\n',
                                      af, dp, altdp, hq, hp, se, bc)

                                df = pd.read_csv(inpt, header=0, index_col=0)
                                df = df[(df.AF >= 0.15) &
                                        (df.Acc.isin(Map.acc))]
                            
                                df["var"] = df.Acc+"_"+\
                                        df.ref+df.pos.astype(str)+df.alt
                                df["type"] = df.apply(lambda x: typer(x),
                                                      axis=1)
                                platform_dict = Map[['acc', 'platform']].\
                                            drop_duplicates().\
                                            set_index('acc').to_dict()
                                df["platform"] = df["Acc"].\
                                                apply(lambda x: 
                                                platform_dict['platform'][x])
                                sample_dict = Map[['acc', 'biosample']].\
                                            drop_duplicates().\
                                            set_index('acc').to_dict()
                                df["biosample"] = df["Acc"].apply(lambda x: 
                                                sample_dict['biosample'][x])
                                df = df[~((df.Group == "NCBI") & (df.DP < 10))]
                                df = df[~(df.DP < dp)]

                                if hq:
                                    df_ = pd.read_csv('TRACE_QA_STATS.csv',
                                                  header=0,
                                                  index_col=None)
                                    df = df[df.Acc.isin(
                                            df_[(df_.pcnt_ref >= 0.5) &
                                                (df_.avg_dp >= 100)].acc)]
                                if hp in [True, False]:
                                    df2_ = pd.read_csv('wuhan_homopolymer.txt',
                                                   sep='\t',
                                                   header=0,
                                                   index_col=None)
                                    spans = [range(int(df2_.at[i, 'start']),
                                               int(df2_.at[i, 'end']))
                                             for i in df2_.index]
                                    pos_ = [p for p in df['pos'].unique()
                                        if any(int(p) in x for x in spans)]
                                    if hp:
                                        df = df[df['pos'].isin(pos_)]
                                    else:
                                        df = df[~df['pos'].isin(pos_)]
                                if se in [True, False]:
                                    if se:
                                        df = df[df.Acc.isin(
                                                    df2[df2.librarylayout ==
                                                                "SINGLE"].acc)]
                                    else:
                                        df = df[(df.platform == 
                                                 "OXFORD_NANOPORE")
                                                | (df.Acc.isin(
                                                    df2[df2.librarylayout !=
                                                            "SINGLE"].acc))]
                                if bc in ['pos', 'bs']:
                                    if bc == 'pos':
                                        df3_ = pd.read_csv('good_match.csv',
                                                      header=0,
                                                      index_col=False)
                                        df3_['Acc'] = ''
                                        df3_.Acc = df3_.biosample.\
                                                apply(lambda x: 
                                                df2[df2.biosample == x].\
                                                acc.to_list())
                                        list_ = df3_[['pos', 'Acc']].\
                                                        values.tolist()
                                        df = df[df.apply(lambda x:
                                                         any((x['pos']==y[0]) &
                                                             (x['Acc'] in y[1])
                                                             for y in list_),
                                                         axis=1)]
                                    if bc == 'bs':
                                        df3_ = pd.read_csv('hq_biosample.csv',
                                                      header=0,
                                                      index_col=False)
                                        df = df[df.Acc.isin(df2[df2.biosample.\
                                                        isin(df3_.biosample)].\
                                                        acc)]
                                df = df[df.AF >= af]
                                df = df[df.DP*df.AF >= altdp]

                                software = {'aligner': {'NCBI': 'HISAT2',
                                                    'BEI': 'BWA-MEM',
                                                    'Galaxy':'BWA-MEM',
                                                    'Gilead': 'SMALT',
                                                    'LANL': 'BWA-MEM',
                                                    'VIR': 'BWA-MEM',
                                                    'Lilly': 'BWA-MEM'},
                                        'caller': {'NCBI': 'GATK',
                                                   'BEI': 'FreeBayes',
                                                   'Galaxy': 'LoFreq',
                                                   'Gilead': 'custom',
                                                   'LANL': 'samtools',
                                                   'VIR': 'LoFreq',
                                                   'Lilly': 'FreeBayes'}}
                                df['Aligner'] = df.Group.apply(lambda x:
                                                        software['aligner'][x])
                                df['Caller'] = df.Group.apply(lambda x:
                                                        software['caller'][x])
    
                                df1 = df
                                df3 = counter(df1, df2, fname_)
                                df4 = counter2(df1)
                                df5 = grouper(df1, df2)
                                df5.to_csv("combined_results"+fname_+".csv")
                                df6 = counter(df5, df2, "vars_"+fname_)
                                scorer(df5, fname_)
                                try:
                                    boxplotter(df5, fname_)
                                except:
                                    pass
                                try:
                                    upsetplotter(df3, fname_)
                                except:
                                    pass
                                try:
                                    upsetplotter(df6, "vars_"+fname_)
                                except:
                                    pass
                                try:
                                    barplotter(df4, fname_)
                                except:
                                    pass
                                try:
                                    scatterplotter(df5, fname_)
                                except:
                                    pass
                                try:
                                    clustermapper(df5, fname_)
                                except:
                                    pass
                                
                            


if __name__ == "__main__":

    import sys

    fname = sys.argv[1]

    df1, df2 = combiner()
    df3 = counter(df1, df2, fname)
    df4 = counter2(df1)
    df5 = grouper(df1, df2)
    df5.to_csv("combined_results"+fname+".csv")
    df6 = counter(df5, df2, "vars_"+fname)
    boxplotter(df5, fname)
    upsetplotter(df3, fname)
    upsetplotter(df6, "vars_"+fname)
    barplotter(df4, fname)
    scatterplotter(df5, fname)
    clustermapper(df5, fname)
