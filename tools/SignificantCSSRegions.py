from quick.webtools.GeneralGuiTool import GeneralGuiTool
from quick.webtools.SplitGtrackIntoLists import splitTrack
from quick.application.ExternalTrackManager import ExternalTrackManager
import time, numpy, sys
from collections import OrderedDict
from quick.util.GenomeInfo import GenomeInfo
from scipy import stats;
#This is a template prototyping GUI that comes together with a corresponding web page.
#

class SignificantCSSRegions(GeneralGuiTool):
    @staticmethod
    def getToolName():
        return "Find significantly diverting regions based on CSS p-values and false discovery rate"

    @staticmethod
    def getInputBoxNames():
        "Returns a list of names for input boxes, implicitly also the number of input boxes to display on page. Each such box will call function getOptionsBoxK, where K is in the range of 1 to the number of boxes"
        #can have two syntaxes: either a list of stings or a list of tuples where each tuple has two items(Name of box on webPage, name of getOptionsBox)
        return [("Reference Genome", "genome"), ("CSS results with p-values", "css"), ("Filtering method", "filter"), ("False discovery rate", "fdr"), ("Number of top scoring regions", "numtop"), ("Maximum distance between significant scores in same region", "wsize")];

    @staticmethod
    def getOptionsBoxGenome():
        return "__genome__"

    @staticmethod    
    def getOptionsBoxCss(prevChoices):
        '''Returns a list of options to be displayed in the first options box
        Alternatively, the following have special meaning:
        '__genome__'
        '__track__'
        ('__history__','bed','wig','...')
        '''
        return ('__history__', 'tabular')

    @staticmethod
    def getOptionsBoxFilter(prevChoices):
        return ["FDR p-value treshold", "top scoring regions"]
        
    @staticmethod
    def getOptionsBoxFdr(prevChoices):
        if prevChoices[2] == "top scoring regions":
            return None;
        return '0.05'

    @staticmethod
    def getOptionsBoxNumtop(prevChoices):
        if prevChoices[2] == "FDR p-value treshold":
            return None;
        return '100'

    @staticmethod
    def getOptionsBoxWsize(prevChoices):
        return '2500'    

    
    #@staticmethod    
    #def getOptionsBox2(prevChoices): 
    #    '''Returns a list of options to be displayed in the second options box, which will be displayed after a selection is made in the first box.
    #    prevChoices is a list of selections made by the web-user in the previous input boxes (that is, list containing only one element for this case)        
    #    '''
    #    return ['']
    
        
    #@staticmethod    
    #def getOptionsBox3(prevChoices):
    #    return '__track__'
        
    #@staticmethod    
    #def getOptionsBox4(prevChoices):
    #    return ['']

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

        genome = choices.genome;
        infile = choices.css;
        fdr_filtering = (choices.filter == "FDR p-value treshold");
        if fdr_filtering:
            FDR = float(choices.fdr)
        else:
            num_top = int(choices.numtop);
        windowSize = int(choices.wsize)

        inFn = ExternalTrackManager.extractFnFromGalaxyTN(infile.split(":"))
        data = open(inFn, "r").read();
        scores, p, addr, windows = cls.preProcessPvalues(data, 2, 3)
        outfile = open(galaxyFn, "w") 
        
        addrs = numpy.array(addr)

        if fdr_filtering:
            # [::-1] --> sorted from smallest to largest
            psorted = numpy.argsort(p)[::-1]
            k = float(len(p))
            n = k
            testp = 0

            #Benjamini-Hochberg procedure
            for pi in psorted:
                if p[pi] <= k/n * FDR:
                    testp = p[pi]
                    break
                k -=1
            # Tuva changed from 1 to 0:
            if k == 0:
                print "NONE FOUND";
                outfile.write("NONE found")
                outfile.close()
                return
        
            print "Pval found:", testp
            filteredaddrs = addrs[p<=testp]
        else:
            scoresorted = numpy.argsort(scores)[::-1];
            scorelimit = scores[scoresorted[num_top-1]];
            filteredaddrs = addrs[scores>=scorelimit];
        
        prevAddr = -10000.
        headers = "##gtrack version: 1.0\n##track type: segments\n##uninterrupted data lines: true\n"+\
                "##no overlapping elements: true\n###seqid\tstart\tend\n"
        outfile.write(headers)
        curchrom = ""
        start = ""
        end = sys.maxint
        prevAddr = -1000000.
        for addr in filteredaddrs:
            addrList = addr.split("\t")
            if addrList[0] != curchrom or int(addrList[1])-windowSize > prevAddr:
                if curchrom != "":
                    newend = prevAddr+windowSize if prevAddr+windowSize < end else end
                    outfile.write(start+"\t"+str(newend)+"\n")
                start = addr
                curchrom = addrList[0]
                end = int(GenomeInfo.getChrLen(genome, curchrom))-1

            prevAddr = int(addr.split("\t")[1])

        newend = prevAddr+windowSize if prevAddr+windowSize < end else end
        outfile.write(start+"\t"+str(newend)+"\n")
        print "Number of regions found", len(filteredaddrs)                                                                                                                 
        if fdr_filtering:
            print "False discoveries", testp*windows
        outfile.close()
    
    
    @classmethod 
    def preProcessPvalues(self, f, scorecol, pcol):
        dicts = []
        addr = []
        lines = f.split("\n")
        scorearray = numpy.zeros(len(lines))+1
        parray = numpy.zeros(len(lines))+1
        for line in lines:
            if (not line) or line[0] == "#":
                continue
            cols = line.split("\t")
            chromaddr = cols[0]+"\t"+cols[1]
            scorearray[len(addr)] = float(cols[scorecol]);
            pval = float(cols[pcol])
            parray[len(addr)] = pval
            addr.append(chromaddr)
        return scorearray[:len(addr)], parray[:len(addr)], addr, len(addr)


    @staticmethod
    def validateAndReturnErrors(choices):
        '''
        Should validate the selected input parameters. If the parameters are not valid,
        an error text explaining the problem should be returned. The GUI then shows this text
        to the user (if not empty) and greys out the execute button (even if the text is empty).
        If all parameters are valid, the method should return None, which enables the execute button.
        '''
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
    #@staticmethod
    #def isDynamic():
    #    return True
    #
    #@staticmethod
    #def getResetBoxes():
    #    return []
    #
    @staticmethod
    def getToolDescription():
        return "<h2>FAQ</h2><h3>What is the aim of this tool?</h3><br />This tool filters out the regions of higher genomic divergence. The input is the Cluster Separation Scores (CSS) from the tool 'Cluster Separation Score'. The regions are filtered with a false discovery rate set by the user, or by selecting the x top scoring CSS regions. The filtering method is selected by the user. "
    
    @staticmethod
    def getFullExampleURL():
        return "https://hyperbrowser.uio.no/comparative/u/tuvakt/h/css-stickleback"
    
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
    #
