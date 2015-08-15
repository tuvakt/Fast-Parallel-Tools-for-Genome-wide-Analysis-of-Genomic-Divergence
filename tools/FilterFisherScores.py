from quick.webtools.GeneralGuiTool import GeneralGuiTool
from quick.application.ExternalTrackManager import ExternalTrackManager
import math, numpy, sys
from quick.util.GenomeInfo import GenomeInfo
from scipy import stats
#This is a template prototyping GUI that comes together with a corresponding web page.
#

class FilterFisherScores(GeneralGuiTool):
    @staticmethod
    def getToolName():
        return "Filter Fisher exact scores into genomic regions"

    @staticmethod
    def getInputBoxNames():
        "Returns a list of names for input boxes, implicitly also the number of input boxes to display on page. Each such box will call function getOptionsBoxK, where K is in the range of 1 to the number of boxes"
        #can have two syntaxes: either a list of stings or a list of tuples where each tuple has two items(Name of box on webPage, name of getOptionsBox)
        return ["Reference genome", "Fisher exact test x% upper quantile results", "Maximum distance between significant scores for a single region", "Quantile of normal distribution to filter", "Percentile of the standard deviation over all windows"]

    @staticmethod
    def getOptionsBox1():
        return "__genome__"
    
    @staticmethod    
    def getOptionsBox2(prevChoices):
        '''Returns a list of options to be displayed in the first options box
        Alternatively, the following have special meaning:
        '__genome__'
        '__track__'
        ('__history__','bed','wig','...')
        '''
        return ('__history__', 'tabular');
    
    @staticmethod    
    def getOptionsBox3(prevChoices): 
        
        '''Returns a list of options to be displayed in the second options box, which will be displayed after a selection is made in the first box.
        prevChoices is a list of selections made by the web-user in the previous input boxes (that is, list containing only one element for this case)        
        '''
        return '100000'

    @staticmethod
    def getOptionsBox4(prevChoices):
        return '0.999'

    @staticmethod
    def getOptionsBox5(prevChoices):
        return '75'

    #@staticmethod
    #def getDemoSelections():
    #    return ['testChoice1','..']
        
    @classmethod    
    def execute(cls, choices, galaxyFn=None, username=''):
        '''Is called when execute-button is pushed by web-user.
        Should print output as HTML to standard out, which will be directed to a results page in Galaxy history.
        If getOutputFormat is anything else than HTML, the output should be written to the file with path galaxyFn.
        If needed, StaticFile can be used to get a path where additional files can be put (e.g. generated image files).
        choices is a list of selections made by web-user in each options box.
        '''
        print 'Executing...'
        genome = choices[0];
        infile = choices[1];
        windowSize = int(choices[2])
        normquantile = float(choices[3])
        percentile = float(choices[4]);
        
        inFn = ExternalTrackManager.extractFnFromGalaxyTN(infile.split(":"));
        data = open(inFn, "r").read()
        fetVals, addr = cls.preProcessPvalues(data, 2)
        stddevs, addr = cls.preProcessPvalues(data, 3)
        output = open(galaxyFn, 'w')
        # Tuva changed sorted elms to FALSE
        output.write("##gtrack version: 1.0\n"+\
                "##track type: segments\n"+\
                "##uninterrupted data lines: true\n"+\
                "##sorted elements: false\n"+\
                "##no overlapping elements: true\n"+\
                "###seqid\tstart\tend\n")
        

        # Calculate limit for FET:
        m = stats.cmedian(fetVals) 
        upperquant = stats.scoreatpercentile(stddevs, percentile)
        qnorm = stats.norm.ppf(normquantile);
        limit = m + qnorm*upperquant
        print "Windows found", sum(fetVals>=limit)
        print 'percentile', percentile, 'normquantile', normquantile
        print "mean",m, 'upperquant', upperquant, 'qnorm', qnorm
        print "Limit", limit 
        addrs = numpy.array(addr)
        filteredaddrs = addrs[fetVals>=limit]
        
        print GenomeInfo.getChrList(genome)

        curchrom = ""
        start = ""
        end = sys.maxint
        prevAddr = -1000000.
        for addr in filteredaddrs:
            addrList = addr.split("\t")
            if addrList[0] != curchrom or int(addrList[1])-windowSize > prevAddr:
                if curchrom != "":
                    newend = prevAddr+windowSize if prevAddr+windowSize < end else end
                    output.write(start+"\t"+str(newend)+"\n")
                start = addr
                curchrom = addrList[0]
                end = int(GenomeInfo.getChrLen(genome, curchrom))-1

            prevAddr = int(addr.split("\t")[1])
        
        newend = prevAddr+windowSize if prevAddr+windowSize < end else end
        output.write(start+"\t"+str(newend)+"\n")
        output.close()
    
    @classmethod 
    def preProcessPvalues(self, f, pcol):
        dicts = []
        addr = []
        lines = f.split("\n")
        parray = numpy.zeros(len(lines))+1
        for line in lines:
            if (not line) or line[0] == "#":
                continue
            cols = line.split("\t")
            chromaddr = cols[0]+"\t"+cols[1]
            pval = float(cols[pcol])
            parray[len(addr)] = pval
            addr.append(chromaddr)
        return parray[:len(addr)], addr
     

    @staticmethod
    def validateAndReturnErrors(choices):
        '''
        Should validate the selected input parameters. If the parameters are not valid,
        an error text explaining the problem should be returned. The GUI then shows this text
        to the user (if not empty) and greys out the execute button (even if the text is empty).
        If all parameters are valid, the method should return None, which enables the execute button.
        '''

        try:
            map(int, choices[2])
        except ValueError:
            return "Only integers allowed for window size"
        try:
            normq = float(choices[4]);
            if normq < 0 or normq > 100:
                raise ValueError;
        except ValueError:
            return "Only values between 0 and 100 allowed for percentile"
        try:
            perc = float(choices[3]);
            if perc < 0 or perc > 1:
                raise ValueError;
        except ValueError:
            return "Only floats between 0 and 1 allowed for qnorm"
        if choices[0] == None:
            return "You must select genome build"
        if choices[0] is not None and (choices[1] == None):
            return "There are no history elements of the right data type, you should probably convert your SNP-data to the correct gtrack format and make sure they are part of your current history."
        return None
    
    @staticmethod
    def isPublic():
        return True
    #
    #@staticmethod
    #def isRedirectTool():
    #    return False
    #
    #@staticmethod
    #def isHistoryTool():
    #    return True
    #
    @staticmethod
    def isDynamic():
        return True
    #
    #@staticmethod
    #def getResetBoxes():
    #    return []
    
    @staticmethod
    def getToolDescription():
        return "<h2>FAQ</h2><h3>What is the aim of this tool?</h3><br />This tool filters out the regions of higher genomic divergence. The input is the Fisher's Exact Test (FET) scores from the tool 'Fisher Exact Test SNP Tool' and the output is a file with the regions with largest genomic divergence. <br><br /> The FET scores are filtered by the following formula, used by <a href='http://www.nature.com/nature/journal/v467/n7315/abs/nature09352.html'>Burke et al. (2010)</a>:<br><br \>median(L10FETx%Q) + qnorm(normquantile)*percentile(stddevs, perc) <br><br> where normquantile (default 0.999) and perc (default 75) is set by the user."
    
    @staticmethod
    def getFullExampleURL():
        return "https://hyperbrowser.uio.no/comparative/u/tuvakt/h/fet-stickleback"
    #@staticmethod
    #def getToolIllustration():
    #    return None
    #
    #@staticmethod
    #def isDebugMode():
    #    return True
    #
    @staticmethod    
    def getOutputFormat(choices):
        '''The format of the history element with the output of the tool.
        Note that html output shows print statements, but that text-based output
        (e.g. bed) only shows text written to the galaxyFn file.
        '''
        return 'gtrack'
    
