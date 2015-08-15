from gold.statistic.MagicStatFactory import MagicStatFactory
from gold.statistic.Statistic import OnlyGloballySplittable, Statistic,MultipleRawDataStatistic, StatisticConcatDictResSplittable, StatisticDynamicDoubleDictSumResSplittable, StatisticConcatResSplittable, StatisticSplittable
from gold.statistic.RawDataStat import RawDataStat
from gold.track.TrackFormat import TrackFormatReq
from quick.track.SlidingWindow import SlidingWindow
from gold.statistic.FormatSpecStat import FormatSpecStat
from gold.util.CustomExceptions import IncompatibleTracksError
from quick.util.FisherCalculator import fisherScorer
from gold.track.GenomeRegion import GenomeRegion
from quick.util.GenomeInfo import GenomeInfo
from gold.result.Results import Results
from collections import OrderedDict
from gold.util.Profiler import Profiler
import numpy, math, sys
from sklearn.manifold import MDS
from scipy import stats
from gold.util.CommonFunctions import smartSum, getClassName, isIter
#from quick.statistic.fisher_cython import fisher_exact_tester;
from quick.statistic.fisher_cython_parallel import fisher_exact_tester;

class FisherExactScoreStat(MagicStatFactory):
    pass



class FisherExactScoreStatUnsplittable(MultipleRawDataStatistic):
    
    def _compute(self):
        reg = self._region;
        
        groupA = self._children[0].getResult()
        groupB = self._children[1].getResult()
        
        # fetching the data as arrays
        apos = groupA.startsAsNumpyArray();
        if len(apos) == 0:
            return {"0":0};
        bpos = groupB.startsAsNumpyArray();
        avals = groupA.valsAsNumpyArray();
        bvals = groupB.valsAsNumpyArray();
        
        regstart = reg.start
        regend = reg.end
        wSize = int(self._kwArgs["wSize"])
        wStep = int(self._kwArgs["wStep"])
        alen = avals.shape[0];
        blen = bvals.shape[0];
       
        # percentile fisher score 
        perc = float(self._kwArgs["percentile"])
        num_win = regend/wStep;
        scores = numpy.zeros(num_win);
        stddev = numpy.zeros(num_win);
 
        print 'region', reg;

        # start the analysis
        fisher_exact_tester(avals, bvals, apos, bpos, regstart, regend, wSize, wStep, alen, blen, perc, scores, stddev);
        print scores;
        print stddev;
        return scores, stddev;

    def _getTrackFormatReq(self):
        return TrackFormatReq(dense=False, allowOverlaps=True)

    #def _createChildren(self):
    #    self._addChild(FormatSpecStat(self._region, self._track, TrackFormatReq(borderHandling=None,allowOverlaps=True,val=None)))
    #    self._addChild(FormatSpecStat(self._region, self._track2, TrackFormatReq(dense=False,val='category', borderHandling=None,allowOverlaps=True )))


