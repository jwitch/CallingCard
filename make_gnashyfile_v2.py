"""
make_gnashyfile.py
written 6/30/16 by RDM
usage
modified 7/8/16 by JNW to look for TTAA in reference
make_gnashyfile -b <barcode filename> -g <genome>  mm or hs default = mm -p <output path>

required
none

not required
-b <barcode filename>, default = ../output_and_analysis/barcodes.txt
-p <path>, default = ../output_and_analysis/  path for input and output

"""

import argparse
import pysam
import pandas as pd
import csv

def sort_gnashy_file(gnashyfilename): ##no changes needed
    gnashy_frame = pd.read_csv(gnashyfilename,delimiter='\t',header=None,names=['chr','pos','reads'])
    gnashy_frame = gnashy_frame.sort_values(['chr','pos'])
    gnashy_frame.to_csv(gnashyfilename,sep='\t',header=False,index=False)
    return [len(gnashy_frame),gnashy_frame.reads.sum()]



##def read_barcode_file(barcode_filename): ##no changes needed
###This function reads in the experiment name, the primer barcode, and the transposon barcodes
###from a file which is in the following format:
###expt name \t primer barcode \t transposon barcode 1,transposon barcode 2, transposon barcode 3 etc.
###It then returns a dictionary with key = a tuple (primer barcode, transposon barcode), value = expt_name
###The last line of the file should not have a return character
##    reader = csv.reader(open(barcode_filename, 'r'),delimiter = '\t')
##    d = {}
##    for row in reader:
##        exname,b1,b2 = row
##        for transposon_barcode in b2.split(","):
##            d[(b1,transposon_barcode)]=exname
##    return d

def make_gnashyfile(filename,outpath,genome):
    #make chromosome list
    if genome == 'hs':
      chr_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
      chr_dict = {'chr1':1,'chr2':2,'chr3':3,'chr4':4,'chr5':5,'chr6':6,'chr7':7,'chr8':8,'chr9':9,'chr10':10,'chr11':11,'chr12':12,'chr13':13,'chr14':14,'chr15':15,'chr16':16,'chr17':17,'chr18':18,'chr19':19,'chr20':20,'chr21':21,'chr22':22,'chrX':23,'chrY':24}
      print("making human gnashyfile")
    elif genome== 'mm':
      chr_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chrX','chrY']
      chr_dict = {'chr1':1,'chr2':2,'chr3':3,'chr4':4,'chr5':5,'chr6':6,'chr7':7,'chr8':8,'chr9':9,'chr10':10,'chr11':11,'chr12':12,'chr13':13,'chr14':14,'chr15':15,'chr16':16,'chr17':17,'chr18':18,'chr19':19,'chrX':20,'chrY':21}
      print("making mouse gnashyfile")
    elif genome == 'ca': ##added Candida
      chr_list = ['Ca21chr1_C_albicans_SC5314','Ca21chr2_C_albicans_SC5314','Ca21chr3_C_albicans_SC5314','Ca21chr4_C_albicans_SC5314','Ca21chr5_C_albicans_SC5314','Ca21chr6_C_albicans_SC5314','Ca21chr7_C_albicans_SC5314','Ca21chrR_C_albicans_SC5314','Ca19-mtDNA']
      chr_dict = {'Ca21chr1_C_albicans_SC5314':1,'Ca21chr2_C_albicans_SC5314':2,'Ca21chr3_C_albicans_SC5314':3,'Ca21chr4_C_albicans_SC5314':4,'Ca21chr5_C_albicans_SC5314':5,'Ca21chr6_C_albicans_SC5314':6,'Ca21chr7_C_albicans_SC5314':7,'Ca21chrR_C_albicans_SC5314':8,'Ca19-mtDNA':9}
      print("making candida gnashyfile")
    #read in experiments and barcodes.  Key = (primer barcode, Transposon barode)
    #Value = expt name
#    barcode_dict = read_barcode_file(bcfilename)
    #initialize quality control dictionary
    qc_dict = {}
    #LOOP THROUGH EXPERIMENTS
    #loop through experiments and make a separate gnashy file for each
#    for expt in list(set(barcode_dict.values())):
      #for each experiment, there will be multiple bam files.  Loop through all of them
      #open output gnashyfile
    expt = filename
    print("Analyzing "+expt)
    output_filename = outpath+expt+".txt"
    output_handle = open(output_filename, 'x') 
      #LOOP THROUGH BAM FILES CORRESPONDING TO 1 experiment
    #for key in barcode_dict.keys(): #this could be made more efficient, but its more clear this way
 #       if barcode_dict[key] == expt:
#          primerBC = key[0]
#          transposonBC = key[1]
#          basename = outpath+expt+"_"+primerBC+"_"+transposonBC
#          sbamFilename = basename+".sorted"
#          pysam.sort(basename+".bam",sbamFilename) ##newest pysam doesn't work; installed pysam==0.8.4 (older versions don't have AlignmentFile attribute)
          #sort and index bamfile
#          sbamFilename = sbamFilename+".bam"
           
#          pysam.index(sbamFilename)
#          print sbamFilename
          #inialize gnashy dictionary
    gnashy_dict = {}
          #make AlignmentFile object
    current_bamfile = pysam.AlignmentFile(filename+".sorted.bam","rb")

          #loop through the chromosomes and pileup start sites
    for chr in chr_list:
        
        aligned_reads_group = current_bamfile.fetch(chr)
            #now loop through each read and pile up start sites

        for aread in aligned_reads_group:
                      #is the read a reverse read?
            if aread.is_reverse:
                      #does it align to a ttaa?
                      #if (aread.query_sequence[-4:]=='TTAA' or aread.query_sequence[-4:]=='ttaa'):
                        #if so, get position and update dictionary
                pos = aread.get_reference_positions()[-1]
                if (chr,pos) in gnashy_dict:
                    gnashy_dict[(chr,pos)]+=1
                else:
                    gnashy_dict[(chr,pos)]=1
            else: #forward read
                      #does it align to a ttaa?
                      #if (aread.query_sequence[0:4]=='TTAA' or aread.query_sequence[0:4]=='ttaa'):
                          #if so, get position and update dicitonary
                pos = aread.get_reference_positions()[0]
                if (chr,pos) in gnashy_dict:
                    gnashy_dict[(chr,pos)]+=1
                else:
                    gnashy_dict[(chr,pos)]=1
                  #output dictionary to gnashy file
    #print(gnashy_dict)
    for key in gnashy_dict:
        output_handle.write("%s\t%s\t%s\n" %(chr_dict[key[0]],key[1], gnashy_dict[key] ))
    output_handle.close()
    #OPEN GNASHY FILE AND SORT BY CHR THEN POS
    qc_dict[expt] = sort_gnashy_file(output_filename)
    #after all experiments have been analyzed, print out qc
    qc_handle = open(outpath+"gnashyQC.txt",'w')
    for key in qc_dict:
      qc_handle.write("%s\t%s\t%s\n" %(key,qc_dict[key][0], qc_dict[key][1] ))
    qc_handle.close()

##if __name__ == '__main__':
##    parser = argparse.ArgumentParser(description='make_gnashyfile_v2.py')
## #   parser.add_argument('-b','--barcodefile',help='barcode filename (full path)',required=False,default='/Users/jessicawitchley/Documents/bowtie2-2.2.9/JW.CCS6_rerun/raw/CC_barcodes.txt')
##    parser.add_argument('-p','--outputpath',help='output path',required=False,default='/Users/jessicawitchley/Documents/bowtie-1.1.2/JW.CCS_10/output_and_analysis')
##    parser.add_argument('-g','--genome',help='genome to be analyzed',required=True)
##    args = parser.parse_args()
##    if not args.outputpath[-1] == "/":
##        args.outputpath = args.outputpath+"/"
##    make_gnashyfile(args.barcodefile,args.outputpath,args.genome)


                            





