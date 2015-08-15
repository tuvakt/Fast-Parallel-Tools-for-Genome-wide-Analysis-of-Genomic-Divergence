# Fast-Parallel-Tools-for-Genome-wide-Analysis-of-Genomic-Divergence
Master's Thesis carried out at the Biomedicial Informatics (BMI) research group at the Department of Informatics at the University of Oslo. 

The thesis can be found at:

In this thesis we have implemented tools for comparative genomics, with two different methods for 
diverting regions in the genome of two populations, a Fisher's Exact Test (FET) and a Cluster Separation Score (CSS). 
For CSS, three different methods for calculating multi-dimensional scaling (MDS) are implemented. The tools are implemented as web tools on the Genomic HyperBrowser, https://hyperbrowser.uio.no/comparative/,
using Python, Cython and C, and the code is parallelized with Pthreads. 

This thesis extends and improves a master's thesis by Torkil Vederhus, that can be found at: 
http://urn.nb.no/URN:NBN:no-39032.

Vederhus (2013) wrote several tools for the Genomic HyperBrowser, (https://hyperbrowser.uio.no/comparative/):
- Fisher Exact Test SNP Tool
- Filter Fisher Scores
- Cluster Separation Score
- Significant CSS Regions
- Convert Stickleback Snps to Gtrack
and two statistics, one for each method, FET and CSS:
- CategoryClusterSeparationStat
- FisherExactScoreStat

and I have made larger and smaller changes to his tools. The largest modificantions are done on the statistics,
while the web tools have smaller changes. 

I have created a new tool, for converting VCF files to the GTrack format used by these tools:
- Convert VCF To Gtrack Tool

All the C and Cython files are created by me. 
