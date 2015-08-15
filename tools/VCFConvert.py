
def convertToGtrackFile(data, pop, genome):

    # output strings
    outdata = "";

    # dictionary of possible 'GT' values
    values = {"./.": -10000, 
              ".|.": -10000,
              "1/0": 0,
              "0/1": 0,
              "1|0": 0,
              "0|1": 0,
              "0/0": 3,
              "0|0": 3,
              "1/1": -3,
              "1|1": -3};

    data = data.strip();
    lines = data.split("\n");    

    # find and process header line
    idx = 0;
    while (not lines[idx].startswith("#CHROM")):
        idx = idx+1;
    header = lines[idx].split("\t");
    chromidx, posidx, formatidx = processHeader(header);
    popidx, pop = getPopulation(header, pop);
   
    # process first line to get index of 'GT' in format
    gtidx = getGT(lines[idx+1].split("\t"), formatidx);
   
    # process lines
    lines = lines[idx+1:];
    for line in lines:
         out = processLine(line, chromidx, posidx, gtidx, popidx, pop, values);
         outdata = outdata + out;

    return outdata;
   
def processHeader(header):
    chromidx = header.index("#CHROM");
    posidx = header.index("POS");
    formatidx = header.index("FORMAT");
    
    return chromidx, posidx, formatidx;


def addHeader(genome):
    return "##gtrack version: 1.0\n"+\
           "##track type: valued points\n"+\
           "##value type: number\n###seqid\t"+\
           ("start\tvalue\tgenomeid\n####genome=%s\n" % genome);
    
def getPopulation(header, pop):
    idx_list = [];
    bad = [];
    for elm in pop:
        try:
            idx_list.append(header.index(elm));
        except:
            print elm, 'not found in file';
            bad.append(elm);

    for elm in bad:
        pop.remove(elm);

    return idx_list, pop;

def getGT(words, formatidx):
    frmat = words[formatidx];
    opts = frmat.split(":");
    return opts.index("GT");

def processLine(line, chromidx, posidx, gtidx, popidx, pop, values):
    words = line.split("\t");
    chrom = words[chromidx];
    pos = words[posidx];
    popvals = [((words[x]).split(":")[gtidx]) for x in popidx];
    out = "";

    for i in range(len(popvals)):
        newline = "%s\t%s\t%s\t%s\n" % (chrom, pos, values[popvals[i]], pop[i]);
        out = out + newline;

    return out;

if __name__ == '__main__':
    # test
    import sys;
    import numpy as np;
    
    if len(sys.argv) < 4:
        print 'Usage: python %s vcf-file pop outfile' % sys.argv[0];
        
    vcf_file = sys.argv[1];
    popfile = sys.argv[2];
    outfilename = sys.argv[3];
    
    pop = [];
    infile = open(popfile, 'r');
    for line in infile:
        pop.append(line.rstrip('\n'));

    print len(pop);
    print pop;
    
    genome = 'test';
    data = open(vcf_file).read();
    
    outfile = open(outfilename, 'w');
    outfile.write(addHeader(genome));
    outfile.write(convertToGtrackFile(data, pop, genome));
    outfile.close();


