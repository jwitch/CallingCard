#usr/bin/python
##########################################################
#Project: Calling_Cards
#Data:2015/03/17
#Author:JINCHUN ZHANG
##########################################################
'''
Used for annotated Significant output
Input: the Significant insertion file, output of the Simpleteste.py
	the refGene data with bed format which can be found at the Code file.
Ouput: the txt file annotated clustered file.
	header:
           Chromosome,
	   Start, 
	   Stop, 
	   BkgHop,
	   BkgHop_withBonus, 
	   ExpHop, 
	   ExpHop_withBonus, 
	   BkgFraction, 
	   ExpFraction, 
	   P_Hyper,
	   P_Poisson,     
           CUTStart,
           CUTEnd,
	   Closest_Upstream_Gene
	   Closest_Upstream_Gene_Common_Name
           CUTstrand,
           CDTStart,
           CDTEnd,
	   Closest_Downstream_Gene
	   Closest_Downstream_Gene_Common_Name	
           CDTstrand
'''


import argparse
import subprocess
import os
from pybedtools import BedTool
import operator
import csv
from genes import RefGenome
from genes import RefGene

parser = argparse.ArgumentParser(prog='Annotated.py', description='Find the closest down-stream and up-stream gene for significant insertions')
parser.add_argument("--a",dest="path", type=str, nargs=1,default=".", help="Type the parent directory of Input and Output")
parser.add_argument("--o",dest="out", type=str, nargs=1,default=".", help="Type the name of Output folder")

args=parser.parse_args()
pwd=args.path[0]
out=args.out[0]
Local=pwd+'Total_Insertions_P_value.txt'
Ref='/Users/jessicawitchley/Documents/bowtie2-2.2.9/JW.CCS_8/Code/C_albicans_genes.corrected.bed'
outfolder=pwd

ref Annotated(pwd,out):
    
    print "Annotated.py Starts!"
    outfolder=pwd+out
    if not os.path.exists(outfolder+'Annotation'):
        os.makedirs(outfolder+'Annotation')

    ##Create file with following headers [Chromosome Start End BkgHop BkgHop_withBonus ExpHop ExpHop_withBonus BkgFraction ExpFraction P_Hyper P_Poisson CUTStart CUTEnd Closest_Upstream_Gene Closest_Upstream_Gene_Common_Name CUTstrand	CDTStart CDTEnd	Closest_Downstream_Gene	Closest_Downstream_Gene_Common_Name CDTstrand]
    with open(outfolder+"Annotation/Total_Insertions_P_value.txt",'r') as testfile:
        test_list=testfile.readlines()
    ##List of lists that will be output to Total_Insertions_Full_Annotation.txt
    formatted_list=[]
    for i in range(1,len(test_list)):
        line=test_list[i]
        line=line.split('\t')
        line[10]=line[10].replace('\n','')
        line[1]=int(line[1]) ##in case sort doesn't like a string of numbers
        formatted_list.append(line)

    sorted_list=formatted_list
    sorted_list.sort(key=lambda row: (row[0],row[1]))
    ##print sorted_list
    with open(outfolder+'/Annotation/Inserts.sorted.bed','w') as inserts:
     ##   inserts=csv.writer(out)
        for line in sorted_list:       
            inserts.writelines(line[0]+'\t'+str(line[1])+'\t'+line[2]+'\n')
    genome=RefGenome(Ref).get_genome()
    ##print genome
    ##annotated_list
    for insert in sorted_list:
        for line in range(1,len(genome)-1): ##does not take into account first and last gene on chromosome
            gene_current=RefGene(genome[line]).get_gene()
            gene_prev=RefGene(genome[line-1]).get_gene()
            gene_next=RefGene(genome[line+1]).get_gene()
            if insert[0] == gene_current.chr:
                ##need to make sure previous gene is on same chromosome
                if insert[0] == gene_prev.chr and insert[0] == gene_next.chr:
                    d_prev=gene_prev.distance_from_start(gene_prev,insert[1])
                    d_current=gene_current.distance_from_start(gene_current,insert[1])
                    d_next=gene_next.distance_from_start(gene_next,insert[1])
                    if d_current < d_prev and d_current < d_next: #find the closest gene
                        if d_prev < d_next: #if next closest gene is upstream report that
                            closest_up_gene=gene_prev
                            d_closest_up=d_prev
                            closest_down_gene=gene_current
                            d_closest_down=d_current
                            break
                        elif d_prev > d_next:
                            closest_up_gene=gene_current
                            d_closest_up=d_current
                            closest_down_gene=gene_next
                            d_closest_down=d_next
                            break
                    elif d_prev < d_current:
                        closest_up_gene=gene_prev
                        d_closest_up=d_prev
                        closest_down_gene=gene_current
                        d_closest_down=d_current
                        print closest_up_gene
                        break
                        
                else:
                    continue
    ##    to_append=str(closest_up_gene.start)+','+str(closest_up_gene.stop)+','+str(closest_up_gene.name)+','+str(closest_up_gene.commonname)+','+str(closest_up_gene.strand)+','+str(d_closest_up)+','+str(closest_down_gene.start)+','+str(closest_down_gene.stop)+','+str(closest_down_gene.name)+','+str(closest_down_gene.commonname)+','+str(closest_down_gene.strand)+','+str(d_closest_down
        ##closest upstream gene start
        insert.append(closest_up_gene.start)
        ##closest upstream gene end
        insert.append(closest_up_gene.stop)
        ##closest upstream gene name
        insert.append(closest_up_gene.name)
        ##closest upstream gene common name
        insert.append(closest_up_gene.commonname)
        ## closest upstream gene strand
        insert.append(closest_up_gene.strand)
        ## closest upstream distance
        insert.append(d_closest_up)
        ##closest downstream gene start
        insert.append(closest_down_gene.start)
        ##closest downstream gene end
        insert.append(closest_down_gene.stop)
        ##closest downstream gene name
        insert.append(closest_down_gene.name)
        ##closest downstream gene common name
        insert.append(closest_down_gene.commonname)
        ## closest downstream gene strand
        insert.append(closest_down_gene.strand)
        ## closest downstream distance
        insert.append(d_closest_down)

    with open(outfolder+'Annotation/Total_Inserts_Full_Annotation.txt','w') as outfile:
        annotation=csv.writer(outfile, delimiter=',')
        annotation.writerow(['Chromosome', 'Start', 'End', 'BkgHop', 'BkgHop_withBonus', 'ExpHop', 'ExpHop_withBonus','BkgFraction',
                             'ExpFraction', 'P_Hyper', 'P_Poisson','CUTStart', 'CUTEnd',
                             'Closest_Upstream_Gene', 'Closest_Upstream_Gene_Common_Name', 'CUTstrand','CUTdistance',
                             'CDTStart', 'CDTEnd','Closest_Downstream_Gene','Closest_Downstream_Gene_Common_Name', 'CDTstrand','CDTdistance'])
        for line in sorted_list:
            annotation.writerow(line)


##inserts=BedTool(outfolder+'/Annotation/Inserts.sorted.bed')
##genes=BedTool(Ref)

##nearest=genes.closest(inserts, d=True,stream=True)
##print nearest

#subprocess.call("sort -k1,1 -k2,2n -k3,3n "+ Local+" | sed '$d' | awk '{print $1,$2,$3}' OFS='\t' > "+outfolder+"/Annotation/Total_Insertions.sorted.bed; echo 'Line1'", shell=True)
#subprocess.call("sort -k1,1 -k2,2n -k3,3n "+ Local+" | sed '$d' > "+outfolder+"/Annotation/Total_Insertions.Sorted.txt; echo 'Line2'", shell=True)
#subprocess.call("closest-features "+outfolder+"/Annotation/Total_Insertions.sorted.bed " +Ref+" > "+outfolder+"/Annotation/Total_Insertions.sorted.annotated.txt; echo 'Line3'" ,shell=True)
#subprocess.call("sed -e 's:NA:NA\tNA\tNA\tNA\tNA\tNA\t:g' -e 's/|/\t/g' "+outfolder+"/Annotation/Total_Insertions.sorted.annotated.txt > "+outfolder+"/Annotation/Total_Insertions.sorted.annotated2.txt; echo 'Line4'", shell=True)
#subprocess.call("cat "+outfolder+"/Annotation/Total_Insertions.sorted.annotated2.txt |awk '{print $5,$6,$7,$8,$9,$11,$12,$13,$14,$15}' OFS='\t' > "+outfolder+"/Annotation/Total_Insertions.sorted.annotated3.txt; echo 'Line5'",shell=True)
#subprocess.call("paste "+outfolder+"/Annotation/Total_Insertions.Sorted.txt "+outfolder+"/Annotation/Total_Insertions.sorted.annotated3.txt > "+outfolder+"/Annotation/Total_Insertions_P_value_Full_Annotated.txt; echo 'Line6'",shell=True )
##subprocess.call("sed -i '1iChromosome	Start	End	BkgHop	BkgHop_withBonus	ExpHop	ExpHop_withBonus	BkgFraction	ExpFraction	P_Hyper	P_Poisson	CUTStart	CUTEnd	Closest_Upstream_Gene	Closest_Upstream_Gene_Common_Name	CUTstrand	CDTStart	CDTEnd	Closest_Downstream_Gene	Closest_Downstream_Gene_Common_Name	CDTstrand' "+outfolder+"/Annotation/Total_Insertions_P_value_Full_Annotated.txt",shell=True )
##subprocess.call("rm "+outfolder+"/Annotation/Total_Insertions.sorted.annotated*.txt "+outfolder+"/Annotation/Total_Insertions.sorted.bed", shell=True)

    print "Annotated.py Completed!"

