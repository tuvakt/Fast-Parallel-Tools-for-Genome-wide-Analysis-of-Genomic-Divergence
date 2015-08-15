from gold.statistic.MagicStatFactory import MagicStatFactory
from gold.statistic.Statistic import OnlyGloballySplittable, Statistic,MultipleRawDataStatistic, StatisticConcatDictResSplittable, StatisticDynamicDoubleDictSumResSplittable, StatisticConcatResSplittable, StatisticSplittable
from gold.statistic.RawDataStat import RawDataStat
from gold.track.TrackFormat import TrackFormatReq
from quick.track.SlidingWindow import SlidingWindow
from gold.statistic.FormatSpecStat import FormatSpecStat
from gold.util.CustomExceptions import IncompatibleTracksError
from quick.webtools.ClusterSeparationScorerStat import separationScore, significanceTreshold, significanceAll
from gold.track.GenomeRegion import GenomeRegion
from quick.util.GenomeInfo import GenomeInfo
from gold.result.Results import Results
from collections import OrderedDict
from gold.util.Profiler import Profiler
import numpy, math
from sklearn.manifold import MDS
from scipy import stats
from gold.util.CommonFunctions import smartSum, getClassName, isIter
#from quick.statistic.css_cython import cluster_separation_scorer;
from quick.statistic.css_cython_parallel import cluster_separation_scorer;

def createSets(k, n):
    if k == 0:
        return [()]
    return ((i,)+s for i in range(n) for s in createSets(k-1, i))

class CategoryClusterSeparationStat(MagicStatFactory):
    pass


class CategoryClusterSeparationStatUnsplittable(MultipleRawDataStatistic):
    
    def _compute(self):
        reg = self._region 
        
        groupA = self._children[0].getResult()
        groupB = self._children[1].getResult()
        
        # fetch data as numpy arrays
        apos = groupA.startsAsNumpyArray();
        if len(apos) == 0:
            return {"0":0};
        bpos = groupB.startsAsNumpyArray();
        avals = groupA.valsAsNumpyArray();
        bvals = groupB.valsAsNumpyArray();
         
        regstart = reg.start;
        regend = reg.end;
        wSize = int(self._kwArgs["wSize"])
        wStep = int(self._kwArgs["wStep"])

        # length of arrays
        alen = avals.shape[0];
        blen = bvals.shape[0]; 

        # monte carlo parameters
        treshold = int(self._kwArgs["mcT"])
        runs = int(self._kwArgs["mcR"])
        
        # choice of distance metric: count differences or average of freq.
        func = self._kwArgs["func"] == "True"
        print 'func', func, type(func);
        if func:
            drosophila = 1;
        else:
            drosophila = 0;
        
        # choice of mds algorithm, 0: cmds, 1: smacof, 2: smacof+cmds
        mds = int(self._kwArgs["mds"]);

        num_win = regend/wStep;
        scores = numpy.zeros(num_win);
        p = numpy.zeros(num_win);
        print '*** region:', reg, '***';
        print 'num_win', num_win;

        # start analysis
        cluster_separation_scorer(avals, bvals, apos, bpos, regstart, regend, wSize, wStep, alen, blen, treshold, runs, drosophila, mds, scores, p);
        print scores;
        print p;
        return scores, p;


    def _getTrackFormatReq(self):
        return TrackFormatReq(dense=False, allowOverlaps=True)

