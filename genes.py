##Script to make class from C. albicans .bed file for comparison to CCS data to find closest gene
##BED file required to be in [Chr,Start,Stop,Gene name,Common name,strand] (Start can be greater than stop if on Crick strand)
##

class RefGenome:
    def __init__(self,bedfile):
        with open(bedfile,'r') as infile:
            genome=infile.readlines()
        bed_list=[]    
        for line in genome:
            line=line.split('\t')
            line[len(line)-1]=line[len(line)-1].replace('\r','') #strip newlines
            line[len(line)-1]=line[len(line)-1].replace('\n','') 
            bed_list.append(line)
            self.genome=bed_list ##list of lists
            
    def get_genome(self):
         return self.genome

class RefGene:
    def __init__(self,gene):
        self.chr = gene[0]
        self.start = gene[1]
        self.stop = gene[2]
        self.name = gene[3]
        self.commonname = gene[4]
        self.strand = gene[5]

    def get_gene(self):
        return self
       
    def distance_from_start(self,gene,binding_event):
        self.distance = abs(int(gene.start)-int(binding_event))
        return self.distance
