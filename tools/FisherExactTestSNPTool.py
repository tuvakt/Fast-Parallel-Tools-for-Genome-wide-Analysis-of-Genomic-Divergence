from quick.webtools.GeneralGuiTool import GeneralGuiTool
import urllib2
from datetime import date
from os import linesep, listdir
import time, re, os, calendar
from gold.util.CommonFunctions import createOrigPath, createDirPath, extractTrackNameFromOrigPath
from quick.util.CommonFunctions import ensurePathExists
from gold.application.GalaxyInterface import GalaxyInterface
from gold.application.StatRunner import AnalysisDefJob
from quick.util.GenomeInfo import GenomeInfo
from gold.track.Track import Track
from gold.track.TrackFormat import TrackFormatReq
from quick.application.ExternalTrackManager import ExternalTrackManager
from gold.origdata.PreProcessTracksJob import PreProcessAllTracksJob
from gold.track.GenomeRegion import GenomeRegion
from quick.util.StaticFile import StaticFile, GalaxyRunSpecificFile


# This is a template prototyping GUI that comes together with a corresponding
# web page.

class FisherExactTestSNPTool(GeneralGuiTool):
    @staticmethod
    def getToolName():
        '''
        Specifies a header of the tool, which is displayed at the top of the
        page.
        '''
        return "Calculate fisher exact test significance score for SNP-positions"

    @staticmethod
    def getInputBoxNames():
        '''
        Specifies a list of headers for the input boxes, and implicitly also the
        number of input boxes to display on the page. The returned list can have
        two syntaxes:
        
            1) A list of strings denoting the headers for the input boxes in
               numerical order.
            2) A list of tuples of strings, where each tuple has
               two items: a header and a key.
        
        The contents of each input box must be defined by the function
        getOptionsBoxK, where K is either a number in the range of 1 to the
        number of boxes (case 1), or the specified key (case 2).
        '''
        return ['Reference genome for species', 'Group 1 SNP data', "Group 2 SNP data", "Sliding Window Size", "Sliding Window Step", "Percentile of L10FET scores in a window", "Output format"] #Alternatively: [ ('box1','1'), ('box2','2') ]

    #@staticmethod
    #def getInputBoxOrder():
    #    '''
    #    Specifies the order in which the input boxes should be displayed, as a
    #    list. The input boxes are specified by index (starting with 1) or by
    #    key. If None, the order of the input boxes is in the order specified by
    #    getInputBoxNames.
    #    '''
    #    return None
    
    @staticmethod    
    def getOptionsBox1(): # Alternatively: getOptionsBoxKey()
        '''
        Defines the type and contents of the input box. User selections are
        returned to the tools in the prevChoices and choices attributes to other
        methods. These are lists of results, one for each input box (in the
        order specified by getInputBoxOrder()).
        
        The input box is defined according to the following syntax:
        
        Selection box:          ['choice1','choice2']
        - Returns: string
        
        Text area:              'textbox' | ('textbox',1) | ('textbox',1,False)
        - Tuple syntax: (contents, height (#lines) = 1, read only flag = False)
        - Returns: string
        
        Password field:         '__password__'
        - Returns: string
        
        Genome selection box:   '__genome__'
        - Returns: string
        
        Track selection box:    '__track__'
        - Requires genome selection box.
        - Returns: colon-separated string denoting track name
        
        History selection box:  ('__history__',) | ('__history__', 'bed', 'wig')
        - Only history items of specified types are shown.
        - Returns: colon-separated string denoting galaxy track name, as
                   specified in ExternalTrackManager.py.
        
        History check box list: ('__multihistory__', ) | ('__multihistory__', 'bed', 'wig')
        - Only history items of specified types are shown.
        - Returns: OrderedDict with galaxy track name as key and selection
                   status (bool) as value.
        
        Hidden field:           ('__hidden__', 'Hidden value')
        - Returns: string
        
        Table:                  [['header1','header2'], ['cell1_1','cell1_2'], ['cell2_1','cell2_2']]
        - Returns: None
        
        Check box list:         OrderedDict([('key1', True), ('key2', False), ('key3', False)])
        - Returns: OrderedDict from key to selection status (bool).
        '''
        return '__genome__'
    
    
    @staticmethod    
    def getOptionsBox2(self): # Alternatively: getOptionsBoxKey()
        return ('__history__', 'gtrack')
        
    @staticmethod    
    def getOptionsBox3(self): # Alternatively: getOptionsBoxKey()
        return ('__history__', 'gtrack')
        
    @staticmethod    
    def getOptionsBox4(self): # Alternatively: getOptionsBoxKey()
        return "2500"

    @staticmethod    
    def getOptionsBox5(self): # Alternatively: getOptionsBoxKey()
        return "500"
    
    @staticmethod
    def getOptionsBox6(self):
        return "0.95"

    @staticmethod
    def getOptionsBox7(self):
        return ["tabular", "html"]
        
    @classmethod    
    def execute(cls, choices, galaxyFn=None, username=''):
        
        from quick.application.GalaxyInterface import GalaxyInterface

        fileformat = choices[6];
        outputFile = open(galaxyFn, "w")

        if fileformat == "html":
            print GalaxyInterface.getHtmlBeginForRuns(galaxyFn)
            print GalaxyInterface.getHtmlForToggles(withRunDescription=False)
            t = calendar.timegm(time.gmtime())
            htmlfile = GalaxyRunSpecificFile(["fet", str(t)], galaxyFn);
            
        genome = choices[0]
        track1 = choices[1].split(":")
        track2 = choices[2].split(":")
        tn1 = ExternalTrackManager.getPreProcessedTrackFromGalaxyTN(genome, track1)
        tn2 = ExternalTrackManager.getPreProcessedTrackFromGalaxyTN(genome, track2)

        windowSize = int(choices[3])
        windowStep = int(choices[4])
        percentile = float(choices[5]);

        #results = {}
        
        # TODO: why this?
        #tr = Track(tn1)
        #tr.addFormatReq(TrackFormatReq(dense=False, allowOverlaps=True))
        
        outputFile.write("#seqid\tstart\tscore\tstddev\n")
        
        if fileformat == "html":
            text = "#seqid\tstart\tscore\tstddev\n";
        print 'chrs:', str(GenomeInfo.getChrList(genome));
        reg = "*"                                                                 
        bins = "*"                                                                                                                                                        
        analysisDef = "Dummy: dummy name ([wStep=%g] [wSize=%g] [percentile=%g])-> FisherExactScoreStat" % (windowStep, windowSize, percentile)                                                       
        userBinSource = GalaxyInterface._getUserBinSource(reg, bins, genome)                                                                                              
        result = GalaxyInterface.runManual([tn1, tn2], analysisDef, reg, bins, genome, galaxyFn=galaxyFn)
        for key in result.getAllRegionKeys():
            chrom = str(key).split(":")[0];
            r = result[key];
            if 'Result' not in r.keys():
                print "skipping chr:", chrom, r;
                continue;
            r = r['Result'];
            scores = r[0];
            stddev = r[1];
            for i in range(len(scores)):
                if (scores[i] != 0):
                    pos = i*windowStep;
                    #if choices[5] == "html":
                        #print "%s\t%s\t%s\t%s\n" % (str(chrom), pos, str(scores[i]), str(stddev[i]))
                    if fileformat == "tabular":
                        outputFile.write("%s\t%s\t%s\t%s\n" % (str(chrom), pos, str(scores[i]), str(stddev[i])))
                    else:
                        text += "%s\t%s\t%s\t%s\n" % (str(chrom), pos, str(scores[i]), str(stddev[i]));
                        
        if fileformat == "html":
            htmlfile.writeTextToFile(text);
            print htmlfile.getLink("Result file");
            print GalaxyInterface.getHtmlEndForRuns()

        outputFile.close()
    
    
    @staticmethod
    def validateAndReturnErrors(choices):
        '''
        Should validate the selected input parameters. If the parameters are not
        valid, an error text explaining the problem should be returned. The GUI
        then shows this text to the user (if not empty) and greys out the
        execute button (even if the text is empty). If all parameters are valid,
        the method should return None, which enables the execute button.
        '''
        try:
            map(int, choices[3:4])
        except ValueError:
            return "Only integers allowed for window size and window step"
        try:
            perc = float(choices[5]);
            if perc < 0 or perc > 1:
                raise ValueError;
        except ValueError:
            return "Only floats between 0 and 1 allowed for percentile"
        if choices[0] == None:
            return "You must select genome build"
        if choices[0] is not None and (choices[1] == None or choices[2] == None):
            return "There are no history elements of the right data type, you should probably convert your SNP-data to the correct gtrack format and make sure they are part\
 of your current history."
        return None 
       
    #@staticmethod
    #def getSubToolClasses():
    #    '''
    #    Specifies a list of classes for subtools of the main tool. These
    #    subtools will be selectable from a selection box at the top of the page.
    #    The input boxes will change according to which subtool is selected.
    #    '''
    #    return None
    #
    @staticmethod
    def isPublic():
    #    '''
    #    Specifies whether the tool is accessible to all users. If False, the
    #    tool is only accessible to a restricted set of users as defined in
    #    LocalOSConfig.py.
    #    '''
        return True
    #
    #@staticmethod
    #def isRedirectTool():
    #    '''
    #    Specifies whether the tool should redirect to an URL when the Execute
    #    button is clicked.
    #    '''
    #    return False
    #
    #@staticmethod
    #def getRedirectURL(choices):
    #    '''
    #    This method is called to return an URL if the isRedirectTool method
    #    returns True.
    #    '''
    #    return ''
    #
    #@staticmethod
    #def isHistoryTool():
    #    '''
    #    Specifies if a History item should be created when the Execute button is
    #    clicked.
    #    '''
    #    return True
    #
    @staticmethod
    def isDynamic():
    #    '''
    #    Specifies whether changing the content of texboxes causes the page to
    #    reload.
    #    '''
        return True
    #
    #@staticmethod
    #def getResetBoxes():
    #    '''
    #    Specifies a list of input boxes which resets the subsequent stored
    #    choices previously made. The input boxes are specified by index
    #    (starting with 1) or by key.
    #    '''
    #    return []
    #
    
    @staticmethod
    def getToolDescription():
    #    '''
    #    Specifies a help text in HTML that is displayed below the tool.
    #    '''
        return '<h2>FAQ</h2><h3>What is the aim of this tool?</h3><br />This tool scores the genomic divergence between two population in sliding windows of the genome based on SNP-data. This is done with the Fishers Exact Test. For each window we calculate a -log10FET score for each position, and select the p% percentile score to represent the window. p is selected by the user. Higher scores means larger divergence. We also calculate the standard deviation of the scores in the window. <br /><br /><h3>How do I format the SNP data for use in this tool?</h3><br /> The headers needed for gtrack files using this tool are:<ul><li>seqid, denoting the chromosome</li><li>start, denoting the position of the SNP</li><li>Value, denoting either the minor allele frequency or the specific individual value for the SNP (3:major, -3: minor, 0:both observed/other, -10000:Not validated)</li><li>genomeid, identifying which individual/population the particular line represents</li></ul>An example of this is:<br />###gtrack version: 1.0<br/>##track type: valued points<br/>##value type: number<br/>###seqid        start   value   genomeid<br/>####genome=gasAcu1<br/>chrVII  104     -3      0<br /><br /><h3> How to only get the scores with the largest genomic divergence? </h3><br />The output can be filtered into significant regions of divergence using the "Filter Fisher Scores"-tool'

    @staticmethod
    def getFullExampleURL():
        return "https://hyperbrowser.uio.no/comparative/u/tuvakt/h/fet-stickleback"

    
    #@staticmethod
    #def getToolIllustration():
    #    '''
    #    Specifies an id used by StaticFile.py to reference an illustration file
    #    on disk. The id is a list of optional directory names followed by a file
    #    name. The base directory is STATIC_PATH as defined by AutoConfig.py. The
    #    full path is created from the base directory followed by the id.
    #    '''
    #    return None
    #
    #@classmethod
    #def isBatchTool(cls):
    #    '''
    #    Specifies if this tool could be run from batch using the batch. The
    #    batch run line can be fetched from the info box at the bottom of the
    #    tool.
    #    '''
    #    return cls.isHistoryTool()
    #
    #@staticmethod
    #def isDebugMode():
    #    '''
    #    Specifies whether debug messages are printed.
    #    '''
    #    return False
    #
    @staticmethod    
    def getOutputFormat(choices):
    #    '''
    #    The format of the history element with the output of the tool. Note
    #    that html output shows print statements, but that text-based output
    #    (e.g. bed) only shows text written to the galaxyFn file.In the latter
    #    case, all all print statements are redirected to the info field of the
    #    history item box.
    #    '''
        format = choices[6]
        if format == "html":
            format = "customhtml"
        return format

