###################################################################
#   SPDI Format parser for added Parsimony between variant callers#
#   Author: Dave Yarmosh, Senior Bioinformatician ATCC 23JUNE2022  #
#                         Version 1.2                             #
###################################################################

import os
import pandas as pd
import re
from Bio import pairwise2

def parsim(fname):
  #handle multiple input data types
  if (fname.split('.')[-1] == 'xlsx') | (fname.split('.')[-1] == '.xls'):
    df = pd.read_excel(fname)
  elif (fname.split('.')[-1] == 'tsv') | (fname.split('.')[-1] == '.txt'):
    df = pd.read_csv(fname,sep='\t')
  elif fname.split('.')[-1] == 'csv':
    df = pd.read_csv(fname)
  else:
    print('Error in input data type')

  #Sort dataframe by group, sample accession, variant position to require that consecutive variants within a group's sample are consecutive
  df = df.sort_values(by=['Group','Acc','pos'])
  fout = fname.split('.')[0].split('/')[-1]+'_parsimonious.csv'
  #set output column names
  #print(df.columns)

  with open(fout, 'w') as f:
    f.write(','.join(df.columns)+'\n')
  #set temporary columns for capturing consecutive indels
  df['group_next'] = df.Group.shift(-1)
  df['Acc_next'] = df.Acc.shift(-1)
  df['pos_next'] = df.pos.shift(-1)
  df['ref_next'] = df.ref.shift(-1)
  df['alt_next'] = df.alt.shift(-1)
  df['type_next'] = df.type.shift(-1)
  df['group_last'] = df.Group.shift(1)
  df['Acc_last'] = df.Acc.shift(1)
  df['pos_last'] = df.pos.shift(1)
  df['ref_last'] = df.ref.shift(1)
  df['alt_last'] = df.alt.shift(1)
  df['type_last'] = df.type.shift(1)
  df['tmp'] = ''
  #make dataframe for deletions.
  altdf = df[df['alt'] == '-']
  #make dataframe for insertions.
  refdf = df[df['ref'] == '-']


  #loop over each row in the SPDI results searching for consecutive deletionss by group and sample accession
  for index,row in altdf.iterrows():
    if ((row.pos +1 == row.pos_next) | (row.pos_last == row.pos -1)) & ((row.Group == row.group_next) | (row.group_last == row.Group)) & ((row.Acc == row.Acc_next) | (row.Acc_last == row.Acc)) & ((row.type_last == 'InDel') | (row.type_next == 'InDel')):
      df.at[index,'tmp'] = True
    else:
      continue

  #loop over each row in the SPDI results searching for consecutive insertions by group and sample accession
  for index,row in refdf.iterrows():
    if ((row.pos +1 == row.pos_next) | (row.pos_last == row.pos -1)) & ((row.Group == row.group_next) | (row.group_last == row.Group)) & ((row.Acc == row.Acc_next) | (row.Acc_last == row.Acc)) & ((row.type_last == 'InDel') | (row.type_next == 'InDel')):
      df.at[index,'tmp'] = True
    else:
      continue

  ##Begin processing consecutive indels
  #collect consecutive variants
  df2 = df[df['tmp'] == True]
  df2 = df2.drop(['tmp'],axis=1)
  #collect nonconsecutive variants
  df3 = df[df['tmp'] != True]
  df3 = df3.drop(['tmp'],axis=1)
 
  ##group consecutive variants by group and sample accession
  df4 = df2.set_index(['Group','Acc','pos','DP','AF','var','type','platform','biosample','Groups','avg_AF','avg_DP','Platforms'])
  #try for the circumstance where none exist
  try:
    #concatenate consecutive variants #consecutive SNPs will be decomposed later
    df4 = df4.groupby(['Group','Acc']).transform(lambda x: ''.join(x)).drop_duplicates().reset_index()
    df4['alt'] = df4['alt'].str.replace(r'-+','-')
    df4['ref'] = df4['ref'].str.replace(r'-+','-')
  except:
    pass
  
  #combine consecutive and nonconsecutive variants back into a single table
  df5 = pd.concat([df3,df4])
  #get rid of temp columns
  df5 = df5.drop(['Acc_next','Acc_last','pos_next','pos_last','group_next','group_last','alt_next','alt_last','ref_next','ref_last'],axis=1)
  #drop duplicates, because the group by will list a record for each consecutive indel
  df5 = df5.sort_values(['Group','Acc','pos','DP'])
  df5 = df5.drop_duplicates(['Group','Acc','pos'],keep='last')
  
  
  ##loop over each row in the SPDI results
  for index,row in df5.iterrows():
      #get unchanging metadata to the left of the position
      leftmeta=','.join([str(index),row['Group'], row['Acc']])
      #get relatively unchanging metadata from the right of the alt allele
      rightmeta=','.join([str(row['DP']),str(row['AF']),row['var'],
                          row['type'],row['platform'],row['biosample'],
                          row['Groups'],str(row['avg_AF']),
                          str(row['avg_DP']),row['Platforms'],row['Aligner'],row['Caller'],'"{0}"'.format(row['Aligners']),'"{0}"'.format(row['Callers'])])
      #nothing needs to be done here. Write to file
      if row['type'] == 'SNP':
         for i in range(0,len(row['ref'])):
           if row['ref'][i] != row['alt'][i]:
             with open(fout, 'a') as f:
               f.write(','.join([leftmeta,
                               str(row['pos']+i),row['ref'][i],row['alt'][i],
                               rightmeta])+'\n')
      #handle well-labeled deletions
      elif row['alt'] == '-':
        with open(fout, 'a') as f:
          f.write(','.join([leftmeta,
                            str(row['pos']),row['ref'],'-',
                            rightmeta])+'\n')
          
      #handle well-labeled insertions
      elif row['ref'] == '-':
        with open(fout, 'a') as f:
          f.write(','.join([leftmeta,
                            str(row['pos']),'-',row['alt'],
                            rightmeta])+'\n')
          
      #remove alt alleles from ref from rightmost position #rightmost preserves the variant position value
      elif ',' in row['alt']: #handle multiple variants at a single site
        for var in row['alt'].split(','): #loop over them
          if len(var) == 1:
            vartype = 'SNP'
          else: #not observed, but handled
            vartype = 'InDel'
            if row['ref'].startswith(var):
  #find latest match between entire alt and any ref alleles
              last_match = row['ref'].find(var) + len(var)
              ref = row['ref'][last_match:] #new
              pos = row['pos'] + last_match #new
              alt = '-'

          
  #remove ref alleles from alt from rightmost position
            elif var.startswith(row['ref']):
              last_match = var.find(row['ref']) + len(row['ref'])
              alt = var[last_match+1:]
              pos = row['pos'] + last_match + 1
              ref = '-'

          rightmeta=','.join([str(row['DP']),str(row['AF']),'"{0}"'.format(row['var']),
                          vartype,row['platform'],row['biosample'],
                          row['Groups'],str(row['avg_AF']),
                          str(row['avg_DP']),row['Platforms'],row['Aligner'],row['Caller'],'"{0}"'.format(row['Aligners']),'"{0}"'.format(row['Callers'])])
          with open (fout,'a') as f:
            f.write(','.join([leftmeta,
                            str(row['pos']),row['ref'],var,rightmeta])+'\n')

      elif row['ref'].startswith(row['alt']):
        #find latest match between entire alt and any ref alleles
        last_match = row['ref'].find(row['alt']) + len(row['alt'])
        ref = row['ref'][last_match:] #new
        pos = row['pos'] + last_match #new
        alt = '-'
        with open(fout, 'a') as f:
          f.write(','.join([leftmeta,
                            str(pos),ref,alt,
                             rightmeta])+'\n')
          
      #remove ref alleles from alt from rightmost position
      elif row['alt'].startswith(row['ref']):
        last_match = row['alt'].find(row['ref']) + len(row['ref'])
        alt = row['alt'][last_match:]
        pos = row['pos'] + last_match
        ref = '-'
        with open(fout, 'a') as f:
          f.write(','.join([leftmeta,
                            str(pos),ref,alt,
                            rightmeta])+'\n')



     #not a SNP, ref does not start with alt, alt does not start with ref, indel is not well-labeled
      else:
          pos=row['pos']
          ref=row['ref']
          alt=row['alt']

          #align alt to ref with match score = 2, mismatch penalty = 0.5, gapopen penalty = 1, gap extension penaly = 0.1 and only save the best alignment
          ms=pairwise2.align.globalms(ref, alt, 2, -.5, -1, -.1,
                                      one_alignment_only=True)

          #loop over the length of the alignment
          for i in range(len(ms[0][0])):
              #if the alleles at position i do not match
              if ms[0][0][i] != ms[0][1][i]:
                  #if there's a deletion
                  if ms[0][0][i] == '-':
                      #Explicitly label InDel
                      rightmeta=','.join([str(row['DP']),str(row['AF']),
                                          row['var'],'InDel',row['platform'],
                                          row['biosample'],row['Groups'],
                                          str(row['avg_AF']),
                                          str(row['avg_DP']),
                                          row['Platforms'],row['Aligner'],row['Caller'],'"{0}"'.format(row['Aligners']),'"{0}"'.format(row['Callers'])])
                  #assume SNP in InDel call #this might be a source of error - no handling for if there's an insertion in this case
                  else:
                      #Explicitly label SNP
                      rightmeta=','.join([str(row['DP']),str(row['AF']),
                                          row['var'],'SNP',row['platform'],
                                          row['biosample'],row['Groups'],
                                          str(row['avg_AF']),
                                          str(row['avg_DP']),
                                          row['Platforms'],row['Aligner'],row['Caller'],'"{0}"'.format(row['Aligners']),'"{0}"'.format(row['Callers'])])
                  ref=ms[0][0][i]
                  alt=ms[0][1][i]
                  with open(fout, 'a') as f:
                    f.write(','.join([leftmeta,str(pos+i),ref,alt,
                                      rightmeta])+'\n')
                  
              else:
                  ref=ms[0][0][i:]
                  #assumes only one set of hyphens/deletions #could put into if for multiple. #Have not witnessed that complicated a variant/alignment
                  alt=ms[0][1][i:].split('-')[0]
                  rightmeta=','.join([str(row['DP']),str(row['AF']),
                                      row['var'],row['type'],row['platform'],
                                      row['biosample'],row['Groups'],
                                      str(row['avg_AF']),
                                      str(row['avg_DP']),
                                      row['Platforms'],row['Aligner'],row['Caller'],'"{0}"'.format(row['Aligners']),'"{0}"'.format(row['Callers'])])
                  if ref.startswith(alt):
                    with open(fout, 'a') as f:
                      f.write(','.join([leftmeta,
                                        str(pos+i+len(alt)),
                                        ref.replace(alt,'',1),
                                        '-',
                                        rightmeta])+'\n')
                      break


if __name__ == "__main__":
    
    import sys
    
    fname = sys.argv[1]
    parsim(fname)
