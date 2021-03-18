#!/usr/bin/python

##########################################################
#Project: Calling_Cards
#Data:2018/01/01
#Author: JNW, modified for python 3
##########################################################

'''
Input: working directory, csv with exp+ref pairs
Wrapper for Calling card seq pipeline
Output: text file of 
*cutoff is used for filter out the independent insertions based on their own count which is the 3rd column of the input dataset
*pwd is the directory that has Reference.txt and Experimental.txt. A21_gnashy_files/Input/
'''
import argparse
import numpy as np
import pandas as pd
import os
import subprocess
import sys
from scipy.stats import poisson
from scipy.stats import hypergeom
from tabulate import tabulate
pd.options.mode.chained_assignment = None
from pybedtools import BedTool
import operator
import csv
from genes import RefGenome
from genes import RefGene
from copy import deepcopy

def CallingCards(pwd,exprefpairs,cutoff,subclustdist,alpha,fake):
    ##creates pandas dataframe for exp+ref pairs
    ##loops through pairs to collect calling card hits
    ##pwd needs to end in /
    ##exprefpairsdf=pd.read_csv(exprefpairs)
    for row in range(0,len(exprefpairs)):
        #print exprefpairs.iloc[row]
        exp=exprefpairs.iloc[row][0]
        ref=exprefpairs.iloc[row][1]
        out=exp+'_'+ref+'_Output/'
        if not os.path.exists(pwd+out):
            os.mkdir(pwd+out)            
        Precutchr(pwd,out,exp,ref,cutoff)
        QTcluster(pwd,out,exp,subclustdist)
        Combine(pwd,out,exp,ref)
        Simpletest(pwd,out,exp,ref,alpha,fake)
        Annotated(pwd,out)

def Precutchr(pwd,out,exp,ref,cutoff):
    print("Precutchr start!")

    for librarygnashy in [str(exp),str(ref)]:
        #print librarygnashy
        precutpath=pwd+out+'PrecutOut/'+librarygnashy
        #print precutpath
        if not os.path.exists(precutpath):
            os.makedirs(precutpath)
        gnashy_infile=pwd+'Input/'+librarygnashy+'.txt'
        #print gnashy_infile
        with open(gnashy_infile,'r') as infile:
                    inf=infile.readlines()
        chr={
            1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[]
            }
        for i in inf:
                #print i
                a=i.split('\t',2)
                #print a
                a[0]=int(a[0]) #chromosome
                a[1]=int(a[1]) #coordinate
                a[2]=int(a[2].replace('\n','')) #num of reads
                if a[2]>cutoff: #remove aberrant reads
                    chr[a[0]].append([a[1],a[2]])
        for i in chr.keys():
            if len(chr[i])>0:
                    chr[i]=np.array(chr[i])
                    uniq=np.unique(chr[i][:,0])
                    new=np.zeros((len(uniq),4),dtype='int32')
                    new[:,1]=uniq
                    new[:,0]=i
                    for n in range(len(uniq)):
                            freq=chr[i][chr[i][:,0]==uniq[n]]
                            freq=freq[:,0]
                            new[n,3]=len(freq)
                            if len(freq)==1:
                                    new[n,2]=1
                            else:
                                    new[n,2]=len(freq)+5**len(freq)	#add bonus for multiple independent tns	
                    outf=precutpath+'/chr'+str(i)+'.txt'
                    with open(outf,'wb') as of:
                            of.write(b'Chromosome\tPosition\tIndependent_Insertions_withBonus\tIndependent_Insertions\n')
                            np.savetxt(of, new,fmt="%d")

    reffilenames = ['chr1.txt', 'chr2.txt', 'chr3.txt', 'chr4.txt', 'chr5.txt', 'chr6.txt', 'chr7.txt', 'chr8.txt', 'chr9.txt']
    with open(pwd+out+'PrecutOut/'+ref+'/combine.txt', 'w+') as outfile:
        ##if chr1.txt file exists...
        if os.path.isfile(pwd+out+'PrecutOut/'+ref+'/chr1.txt'):
            with open(pwd+out+'PrecutOut/'+ref+'/chr1.txt', 'r') as infile:
                #print infile.readline().strip()
                outfile.write(infile.readline().strip())
            for fname in reffilenames:
                if os.path.exists(pwd+out+'PrecutOut/'+ref+'/'+fname):
                    with open(pwd+out+'PrecutOut/'+ref+'/'+fname, 'r') as infile:
                        for line in infile.readlines()[1:]:
                            #print line
                            outfile.write(line)
        ##if chr2.txt file exists...
        elif os.path.isfile(pwd+out+'PrecutOut/'+ref+'/chr2.txt'):
            with open(pwd+out+'PrecutOut/'+ref+'/chr2.txt', 'r') as infile:
                #print infile.readline().strip()
                outfile.write(infile.readline().strip())
            for fname in reffilenames:
                if os.path.exists(pwd+out+'PrecutOut/'+ref+'/'+fname):
                    with open(pwd+out+'PrecutOut/'+ref+'/'+fname, 'r') as infile:
                        for line in infile.readlines()[1:]:
                            #print line
                            outfile.write(line)
        ##if chr3.txt file exists...
        elif os.path.isfile(pwd+out+'PrecutOut/'+ref+'/chr3.txt'):
            with open(pwd+out+'PrecutOut/'+ref+'/chr3.txt', 'r') as infile:
                #print infile.readline().strip()
                outfile.write(infile.readline().strip())
            for fname in reffilenames:
                if os.path.exists(pwd+out+'PrecutOut/'+ref+'/'+fname):
                    with open(pwd+out+'PrecutOut/'+ref+'/'+fname, 'r') as infile:
                        for line in infile.readlines()[1:]:
                            #print line
                            outfile.write(line)
        ##if chr4.txt file exists...
        elif os.path.isfile(pwd+out+'PrecutOut/'+ref+'/chr4.txt'):
            with open(pwd+out+'PrecutOut/'+ref+'/chr4.txt', 'r') as infile:
                #print infile.readline().strip()
                outfile.write(infile.readline().strip())
            for fname in reffilenames:
                if os.path.exists(pwd+out+'PrecutOut/'+ref+'/'+fname):
                    with open(pwd+out+'PrecutOut/'+ref+'/'+fname, 'r') as infile:
                        for line in infile.readlines()[1:]:
                            #print line
                            outfile.write(line)
        ##if chr5.txt file exists...
        elif os.path.isfile(pwd+out+'PrecutOut/'+ref+'/chr5.txt'):
            with open(pwd+out+'PrecutOut/'+ref+'/chr5.txt', 'r') as infile:
                #print infile.readline().strip()
                outfile.write(infile.readline().strip())
            for fname in reffilenames:
                if os.path.exists(pwd+out+'PrecutOut/'+ref+'/'+fname):
                    with open(pwd+out+'PrecutOut/'+ref+'/'+fname, 'r') as infile:
                        for line in infile.readlines()[1:]:
                            #print line
                            outfile.write(line)
        ##if chr6.txt file exists...
        elif os.path.isfile(pwd+out+'PrecutOut/'+ref+'/chr6.txt'):
            with open(pwd+out+'PrecutOut/'+ref+'/chr6.txt', 'r') as infile:
                #print infile.readline().strip()
                outfile.write(infile.readline().strip())
            for fname in reffilenames:
                if os.path.exists(pwd+out+'PrecutOut/'+ref+'/'+fname):
                    with open(pwd+out+'PrecutOut/'+ref+'/'+fname, 'r') as infile:
                        for line in infile.readlines()[1:]:
                            #print line
                            outfile.write(line)
        ##if chr7.txt file exists...
        elif os.path.isfile(pwd+out+'PrecutOut/'+ref+'/chr7.txt'):
            with open(pwd+out+'PrecutOut/'+ref+'/chr7.txt', 'r') as infile:
                #print infile.readline().strip()
                outfile.write(infile.readline().strip())
            for fname in reffilenames:
                if os.path.exists(pwd+out+'PrecutOut/'+ref+'/'+fname):
                    with open(pwd+out+'PrecutOut/'+ref+'/'+fname, 'r') as infile:
                        for line in infile.readlines()[1:]:
                            #print line
                            outfile.write(line)
        ##if chr8.txt file exists...
        elif os.path.isfile(pwd+out+'PrecutOut/'+ref+'/chr8.txt'):
            with open(pwd+out+'PrecutOut/'+ref+'/chr8.txt', 'r') as infile:
                #print infile.readline().strip()
                outfile.write(infile.readline().strip())
            for fname in reffilenames:
                if os.path.exists(pwd+out+'PrecutOut/'+ref+'/'+fname):
                    with open(pwd+out+'PrecutOut/'+ref+'/'+fname, 'r') as infile:
                        for line in infile.readlines()[1:]:
                            #print line
                            outfile.write(line)
        
'''Define Function _SUBCLUST_:
	INPUT:  
        index:
		pos:
		d: universal variate, the max diameter
	OUTPUT: subclusters of input  points and  within distance le 2.5 kb
'''
def _SUBCLUST_ (index,pos,subclustdist):
	aa=0
	pos=sorted(pos)
	while aa==0:
		left_index=min(index)
		right_index=max(index)
		sub=[pos[x] for x in index] 
		if left_index==0 and right_index<len(pos)-1 and pos[right_index+1]-pos[left_index]<=subclustdist:
			sub.append(pos[right_index+1])	
			index.append(right_index+1)
		elif right_index==len(pos)-1 and left_index>0 and pos[right_index]-pos[left_index-1]<=subclustdist:
			index.append(left_index-1)
			sub.append(pos[left_index-1])
		elif left_index>0 and right_index<len(pos)-1:
			left=pos[left_index-1] 
			right=pos[right_index+1]
			center=(pos[right_index]+pos[left_index])/2
			if center-left>right-center and right-pos[left_index]<=subclustdist:
                                sub.append(right)
                                index.append(right_index+1)
			elif center-left<right-center and pos[right_index]-left<=subclustdist:
                                sub.append(left)
                                index.append(left_index-1)
			elif center-left==right-center and right-left<=subclustdist:
                                sub.append(right)
                                sub.append(left)
                                index.append(left_index-1)
                                index.append(right_index+1)
			else: aa=1
		else: aa=2
	return sub

'''Define Function _MAX_:
	INPUT: position list
	OUTPUT: one subcluster with maximun within distance
''' 

def _MAX_ (pos,subclustdist):
	MAX=0
	MAX_sub=[]
	for p in range(len(pos)):
		index=[p]
		subp=_SUBCLUST_(index,pos,subclustdist)
		if max(subp)-min(subp)>=MAX:
                    MAX=max(subp)-min(subp)
                    MAX_sub=subp
	return MAX_sub

'''Define Function _QTCLUST_
	INPUT: pos:position set
	OUTPUT: subclustered data with format 'position number_of_cluster'
		ps. 99999 refers to singleton 
'''
def _QTCLUST_(pos,subclustdist):
        subnum=0
        pp=0
        qtclustered=[]
        while pp==0:
                pos=sorted(pos)
                length=len(pos)
                if length==0: 
                        pp=1
                if length==1:
                        pp=2
                        qtclustered.append([pos[0],99999])
                elif length>1:
                        sub=_MAX_(pos,subclustdist)
                        if len(sub)==1:
                                qtclustered.append([sub[0],99999])
                        else:
                                subnum=subnum+1 
                                for b in sub:
                                        qtclustered.append([b,subnum])
                        pos=list(set(pos)-set(sub))
                        pos=sorted(pos)
        return qtclustered

'''Read one Chr once and process
	Input: Chr sorted by pos with format 'chr position count'
Output: save QTclustered dataframe data for each chr with format'Chr,Position,Window,Count'
'''

def QTcluster(pwd,out,exp,subclustdist):
        print("QTcluster.py start!")
        print("The maximum distance for a subcluster is: ", subclustdist)

        if not os.path.exists(pwd+out+'PrecutOut/'+exp):
                print("Wrong Directory, please put the one that contains the experiment data")
                sys.exit("Ending Script")	    

        if not os.path.exists(pwd+out+'QTout'):
                os.makedirs(pwd+out+'QTout')
        for i in range(1,9):
                inf=str(pwd)+str(out)+'PrecutOut/'+str(exp)+'/chr'+str(i)+'.txt'
                if os.path.isfile(inf):
                        chr=np.loadtxt(inf,skiprows=1,dtype=int)
                        #sfprint 'chr'+str(i)
                        if len(chr.shape)==2:
                                positions=chr[:,1]
                                qtc=np.array(_QTCLUST_(positions,subclustdist))
                                df=pd.DataFrame(qtc,index=qtc[:,0],columns=['Position','Cluster'])
                                s1=pd.Series(chr[:,2],index=chr[:,1])
                                s2=pd.Series(chr[:,3],index=chr[:,1])
                                df['Independent_Insertion_withBonus']=s1
                                df['Independent_Insertion']=s2
                                df['Chromosome']=i
                                df=df.sort_values(by='Position',axis=0)	
                        else:
                                df=pd.DataFrame(columns=["Position","Cluster","Independent_Insertion_withBonus","Independent_Insertion","Chromosome"])
                                df['Chromosome']=[chr[0]]
                                df["Position"]=[chr[1]]
                                df["Cluster"]=['99999']
                                df["Independent_Insertion_withBonus"]=[chr[2]]
                                df["Independent_Insertion"]=[chr[3]]
                        outf=pwd+out+'QTout/chr'+str(i)+'clustered.csv'
                        df.to_csv(outf, index=False,header="Position,Cluster,Independent_Insertion_withBonus,Independent_Insertion,Chromosome") 

        print("QTcluster completed!")

def Combine(pwd,out,exp,ref):
        print("Combine.py started")
        outpath=pwd+out
        if not os.path.exists(outpath):
                print("Wrong Directory, please put the one ended with Output/")
                sys.exit("Ending Script")
        outf=outpath+'QTout/'+ref+'_Combine.csv'
        #print outf
        counter=0
        for chromosome in range(1,9): ##changed for C albicans
                expf=outpath+'QTout/chr'+str(chromosome)+'clustered.csv'
                #print expf
                bkgf=outpath+'PrecutOut/'+ref+'/chr'+str(chromosome)+'.txt'
                #print bkgf
                if os.path.isfile(expf):
                        #print "True"
                        counter=counter+1
                        exp=pd.read_csv(expf)
                        clust=np.array(exp['Cluster']).tolist()
                        bkg=np.loadtxt(bkgf,skiprows=1,dtype=int)
                        comb=_COMBINE_(exp,bkg,chromosome)
                        if len(np.unique(clust))==1:
                                comb=_TRANS1_(comb,chromosome)
                        if len(np.unique(clust))>1:
                                clust=np.unique(clust)[-2]
                                comb=_FILL_(comb,clust)
                                comb=_TRANS_(comb,clust,chromosome)
                        comb=comb.astype('int64')
                        if counter==1:
                                comb.to_csv(outf,header=True,index=False)
                        else:
                                comb.to_csv(outf,header=False,index=False,mode='a')
                        counter=counter+1
        print('Combine.py Completed!')



def _COMBINE_ (exp,bkg,chromosome):
        exppos=exp['Position'].tolist()
        window=exp['Cluster'].tolist()
        window=pd.Series(window,index=exppos)
        expcount=exp['Independent_Insertion_withBonus'].tolist()
        expcount2=exp['Independent_Insertion'].tolist()
        expcount=pd.Series(expcount,index=exppos)
        expcount2=pd.Series(expcount2,index=exppos)
        bkgcount=pd.Series(list(bkg[:,2]),index=bkg[:,1])
        bkgcount2=pd.Series(list(bkg[:,3]),index=bkg[:,1])
        pos=np.append(np.array(exppos),bkg[:,1])
        pos=np.unique(pos)
        pos=pd.Series(list(pos),index=list(pos))
        new=pd.DataFrame({'Position':pos,'Cluster':window,'ExpHop_withBonus':expcount,'BkgHop_withBonus':bkgcount,'ExpHop':expcount2,'BkgHop':bkgcount2})
        new.index=range(len(new))
        return new

def _FILL_ (comb,clust):
        for p in range(1,clust+1):
                ind=comb[comb['Cluster']==p].index
                if len(ind)>2 and (ind[-1]-ind[0]-len(ind)) != -1:
                        for a in range(ind[0],ind[-1]+1):
                                comb.loc[a,'Cluster']=p
                comb=comb.fillna(0)
        return 	comb

def _TRANS_ (comb,clust,chromosome):
        data=comb[comb['Cluster']==99999]
        data.insert(0,'Chromosome',chromosome)
        s=data['Cluster']
        del data['Cluster']
        data.insert(1,'Cluster',s)
        data.insert(2,'Start',data['Position'])
        data.insert(3,'Stop',data['Position'])
        del data['Position'] 
        if not data.empty:
            idx=data.index[-1]
            data=data.astype('int64')
            for grp in range(1,clust+1):
                    block=comb[comb['Cluster']==grp]
                    N=np.array(block['Position'])
                    data.loc[idx+grp,'Chromosome']=chromosome
                    data.loc[idx+grp,'Start']=np.amin(N)
                    data.loc[idx+grp,'Stop']=np.amax(N)
                    data.loc[idx+grp,'BkgHop']=sum(block['BkgHop'])
                    data.loc[idx+grp,'ExpHop']=sum(block['ExpHop'])
                    data.loc[idx+grp,'BkgHop_withBonus']=sum(block['BkgHop_withBonus'])
                    data.loc[idx+grp,'ExpHop_withBonus']=sum(block['ExpHop_withBonus'])
                    data.loc[idx+grp,'Cluster']=grp
        data.index=range(len(data))
        return data

def _TRANS1_(comb,chromosome):
        comb=comb.fillna(0)
        data=comb[comb['Cluster']==99999]
        data.insert(0,'Chromosome',chromosome)
        s=data['Cluster']
        del data['Cluster']
        data.insert(1,'Cluster',s)
        data.insert(2,'Start',data['Position'])
        data.insert(3,'Stop',data['Position'])
        del data['Position']
        if not data.empty:
            idx=data.index[-1]
        data=data.astype('int64')
        return data


def Simpletest(pwd,out,exp,ref,alpha,fake):
        print("Simpletest started!")
        
        if not os.path.exists(pwd+out+'QTout'):
                print("Wrong Directory, please put the one contains a folder QTout which has a file named Combine.csv.")
                sys.exit("Ending Script")
        if not os.path.exists(pwd+out+'Annotation'):
                os.makedirs(pwd+out+'Annotation')
        summary=[]
        outfile_name=pwd+out+'Annotation/Total_Insertions_P_value.txt'
        infile_name=pwd+out+'QTout/'+ref+'_Combine.csv'
        #print "infile_name"
        bkgtotal=pwd+out+'PrecutOut/'+ref+'/combine.txt'
        bkgtotal=np.loadtxt(bkgtotal,dtype=int,skiprows=1)
        #print bkgtotal
        bkgtotal=sum(bkgtotal[:,3])
        ALL=pd.read_csv(infile_name)
        exptotal=sum(ALL['ExpHop'])
        counter=0

        for c in range(1,9):
                #print c
                block=ALL[ALL['Chromosome']==c]
                if len(block)>0:
                        counter=counter+1
                        #print counter
                        out=_TEST_(block,exptotal,bkgtotal,alpha,float(fake))	
                        block=out[0]
                        x='chr'+str(c)
                        block['Chromosome']=x
                        singleton=block[block['Cluster']==99999]
                        singleton=len(singleton[singleton['ExpHop']==1])
                        summary.append([x,(len(block)-singleton),singleton,out[1],out[2]])
                        del block['Cluster']
                        if counter==1:
                            block.to_csv(outfile_name,sep='\t',header=True,index=False)
                        else:
                            #print outfile_name
                            block.to_csv(outfile_name,sep='\t',header=False,index=False,mode='a')
        total=['Total',0,0,0,0]
        for item in summary:
                total[1]=total[1]+item[1]
                total[2]=total[2]+item[2]
                total[3]=total[3]+item[3]
                total[4]=total[4]+item[4]

        summary.append(total)
        print("Simpletest Result Summary:")
        print("Reference total Hops:", bkgtotal)
        print("Experiment Total Hops:", exptotal)
        print(tabulate(summary,headers=['Chr','Cluster','Singleton','Sigificant_By_Poisson_CDF','Significant_By_Hypergeometric_CDF'],tablefmt='orgtbl'))
        print("Simpletest.py Completed!")


'''
   Function _TEST_: Test the Null Hypothesis that the SP1 == WT.
	#1.Poisson distribution.
	# Test statistics: Pr(X>=SP1hop|WT)
	#2.Hypergeometirc districution
        # scistat.hypergeom.cdf(x,M,n,N)
        # where x is observed number of type I events (white balls in draw) (experiment hops at locus)
        # M is total number of balls (total number of hops,wt+exp total)
        # n is total number of white balls (total number of experimental hops)
        # N is the number of balls drawn (total hops at a locus)

'''
def _TEST_ (block,exptotal,bkgtotal,alpha,pseudo_count):
        sigPoi=0
        sigHyp=0
        ind=block.index
        for i in block.index:
                lamda=block.loc[i,'BkgHop_withBonus']
                obs=block.loc[i,'ExpHop_withBonus']
                BkgHop=block.loc[i,'BkgHop']
                ExpHop=block.loc[i,'ExpHop']
                P_poisson=1-poisson.cdf(int(obs)-1,lamda+pseudo_count)
                if int(BkgHop) < 100000000 and int(ExpHop) < 100000000:
                        P_hyper=1-hypergeom.cdf(ExpHop-1,(bkgtotal+exptotal),exptotal, (ExpHop+BkgHop))
                else:
                        P_hyper='***'
                block.loc[i,'BkgFraction']=float(BkgHop/bkgtotal)
                block.loc[i,'ExpFraction']=float(ExpHop/exptotal)
                block.loc[i,'P_Hyper']= P_hyper
                block.loc[i,'P_Poisson']=P_poisson
                if P_poisson < alpha and lamda*obs != 0:
                        sigPoi=sigPoi+1
                if P_hyper < alpha:
                        sigHyp=sigHyp+1
        return block, sigPoi, sigHyp

def Annotated(pwd,out):
    print("Annotated.py Starts!")
    outfolder=pwd+out
    if not os.path.exists(outfolder+'Annotation'):
        os.makedirs(outfolder+'Annotation')
    Ref=pwd+'genome/C_albicans_genes.bed'
    ##Create file with following headers [Chromosome Start End BkgHop BkgHop_withBonus ExpHop ExpHop_withBonus BkgFraction ExpFraction P_Hyper P_Poisson CUTStart CUTEnd Closest_Upstream_Gene Closest_Upstream_Gene_Common_Name CUTstrand	CDTStart CDTEnd	Closest_Downstream_Gene	Closest_Downstream_Gene_Common_Name CDTstrand]
    with open(outfolder+"Annotation/Total_Insertions_P_value.txt",'r') as testfile:
        test_list=testfile.readlines()
    ##print test_list
    ##List of lists that will be output to Total_Insertions_Full_Annotation.txt
    formatted_list=[]
    for i in range(1,len(test_list)):
        line=test_list[i]
        line=line.split('\t')
        line[10]=line[10].replace('\r','') ##double check last element
        line[10]=line[10].replace('\n','')
        line[1]=int(line[1]) ##in case sort doesn't like a string of numbers
        formatted_list.append(line)

    sorted_list=formatted_list
    #print(sorted_list)
    sorted_list.sort(key=lambda row: (row[0],row[1]))
    ##print sorted_list
    with open(outfolder+'/Annotation/Total_Insertions.sorted.bed','w') as inserts:
        for line in sorted_list:       
            inserts.writelines(line[0]+'\t'+str(line[1])+'\t'+line[2]+'\n')

    with open(outfolder+'/Annotation/Total_Insertions.sorted.txt','w') as inserts:
        for line in sorted_list:
            inserts.writelines(line[0]+'\t'+str(line[1])+'\t'+line[2]+'\t'+line[3]+'\t'+line[4]+'\t'+line[5]+'\t'+line[6]+'\t'+line[7]+'\t'+line[8]+'\t'+line[9]+'\t'+line[10]+'\n')

    subprocess.call("closest-features --dist "+outfolder+"/Annotation/Total_Insertions.sorted.bed " +Ref+" > "+outfolder+"/Annotation/Total_Insertions.sorted.annotated.txt" ,shell=True)

    f = open(outfolder+"/Annotation/Total_Insertions.sorted.annotated.txt",'r')
    filedata = f.read()
    f.close()

    newdata = filedata.replace("\n|","\t")
    newdata = newdata.replace("|","\t")

    f = open(outfolder+"/Annotation/Total_Insertions.sorted.annotated.replace.txt",'w')
    f.write(newdata)
    f.close()
    
    subprocess.call("sed -e 's:NA:NA\tNA\tNA\tNA\tNA\tNA\t:g' -e 's/|/\t/g' "+outfolder+"/Annotation/Total_Insertions.sorted.annotated.replace.txt > "+outfolder+"/Annotation/Total_Insertions.sorted.annotated2.txt", shell=True)
    subprocess.call("cat "+outfolder+"/Annotation/Total_Insertions.sorted.annotated2.txt |awk '{print $5,$6,$7,$8,$9,$10,$12,$13,$14,$15,$16,$17}' OFS='\t' > "+outfolder+"/Annotation/Total_Insertions.sorted.annotated3.txt",shell=True)
    subprocess.call("paste "+outfolder+"/Annotation/Total_Insertions.Sorted.txt "+outfolder+"/Annotation/Total_Insertions.sorted.annotated3.txt > "+outfolder+"/Annotation/Total_Insertions_P_value_Full_Annotation.txt",shell=True )
#    subprocess.call("sed -i '1iChromosome	Start	End	BkgHop	BkgHop_withBonus	ExpHop	ExpHop_withBonus	BkgFraction	ExpFraction	P_Hyper	P_Poisson	CUTStart	CUTEnd	Closest_Upstream_Gene	Closest_Upstream_Gene_Common_Name	CUTstrand	CDTStart	CDTEnd	Closest_Downstream_Gene	Closest_Downstream_Gene_Common_Name	CDTstrand' "+outfolder+"/Annotation/Total_Insertions_P_value_Full_Annotation.txt",shell=True )
 #   subprocess.call("rm "+outfolder+"/Annotation/Total_Insertions.sorted.annotated*.txt "+outfolder+"/Annotation/Total_Insertions.sorted.bed", shell=True)

##    genome=RefGenome(Ref).get_genome()

##    intergenics=[] ##create list of intergenic regions
##    for idx in range(0,len(genome)-1):
##        #print(idx)
##        intergenics.append(genome[idx][:]+genome[idx+1][:])
##        
##    cluster_assignment=deepcopy(sorted_list)
##    for insert in cluster_assignment:
##        for line in intergenics: #now takes into account ends of chromosomes
##                if insert[0] == line[0] and line[0] == line[6]: #make sure data on same chromosome
##                    insertStart=int(insert[1])
##                    insertEnd=int(insert[2])
##                    intergenic5start=int(line[1])
##                    intergenic5end=int(line[2])
##                    intergenic3start=int(line[7])
##                    intergenic3end=int(line[8])
##                    if (insertStart > intergenic5start and insertEnd < intergenic3end): #make sure data is bounded
##                        insert.extend(line[1:6])
##                        insert.append(insertStart-intergenic5end)
##                        insert.extend(line[7:12])
##                        insert.append(intergenic3start-insertStart)
##                        ##can result in two different intergenic regions concatenated if insert is in the middle of a gene
##
##                    #print(insert)
##    cluster_assignment_cut=[]
##    for i in range(0,len(cluster_assignment)):
##        insert=cluster_assignment[i]
##        #print(insert)
##        if len(insert) > 23 and len(insert) < 36:
##            #print(insert)
##            if abs(int(insert[16])) < abs(int(insert[34])):
##                firstinsert=insert[0:23]
##                cluster_assignment_cut.append(firstinsert)
##            elif abs(int(insert[16])) > abs(int(insert[34])):
##                #print(insert)
##                secondinsert=insert[0:11]
##                secondinsert.extend(insert[23:])
##                #print(secondinsert)
##                cluster_assignment_cut.append(secondinsert)
##        #elif len(insert) > 29 and len(insert) < 40:
##            #print(insert)
##        else:
##            cluster_assignment_cut.append(insert)

##    filenohead=pd.read_csv(outfolder+"/Annotation/Total_Insertions_P_value_Full_Annotation.txt",sep="\t",header=None)
##    filehead=filenohead.to_csv(outfolder+"/Annotation/Total_Insertions_P_value_Full_Annotation.txt",sep="\t",header=['Chromosome', 'Start', 'End', 'BkgHop', 'BkgHop_withBonus', 'ExpHop', 'ExpHop_withBonus','BkgFraction',
##                             'ExpFraction', 'P_Hyper', 'P_Poisson','CUTStart', 'CUTEnd',
##                             'Closest_Upstream_Gene', 'Closest_Upstream_Gene_Common_Name', 'CUTstrand','CUTdistance',
##                             'CDTStart', 'CDTEnd','Closest_Downstream_Gene','Closest_Downstream_Gene_Common_Name', 'CDTstrand','CDTdistance'])
##    with open(outfolder+'Annotation/Total_Inserts_Full_Annotation.txt','w') as outfile:
##        annotation=csv.writer(outfile, delimiter=',')
##        annotation.writerow(['Chromosome', 'Start', 'End', 'BkgHop', 'BkgHop_withBonus', 'ExpHop', 'ExpHop_withBonus','BkgFraction',
##                             'ExpFraction', 'P_Hyper', 'P_Poisson','CUTStart', 'CUTEnd',
##                             'Closest_Upstream_Gene', 'Closest_Upstream_Gene_Common_Name', 'CUTstrand','CUTdistance',
##                             'CDTStart', 'CDTEnd','Closest_Downstream_Gene','Closest_Downstream_Gene_Common_Name', 'CDTstrand','CDTdistance'])
####        for line in cluster_assignment_cut:
##            annotation.writerow(line)
            
    print("Annotated.py Completed!")
