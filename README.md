# Fast-Parallel-Tools-for-Genome-wide-Analysis-of-Genomic-Divergence
Master's Thesis carried out at the Biomedicial Informatics (BMI) research group at the Department of Informatics at the University of Oslo. 

Comparative genomics is useful for finding the evolutionary relationships between different organisms. 
Regions of parallel evolution can be identified, and we can gain a better understanding of evolution and the development of different species. 
Tools for comparative genomics have previously been developed by \citeauthor{torkil} (\citeyear{torkil}).
The tools, however, had some shortcomings, such as the usability and the long runtime of the tools.
The aim of this thesis was to improve the tools for comparative genomics.

We have implemented two different methods for locating diverting regions in the genome of two populations, a Fisher's Exact Test (FET) and a Cluster Separation Score (CSS). 
For CSS, three different methods for calculating multi-dimensional scaling (MDS) were implemented. The tools were implemented as web tools on the Genomic HyperBrowser, using Python, Cython and C,
and the code was parallelized with Pthreads. 
We have simplified the file format and given the user more choices in parameters, and created a tool for converting VCF files to our file format.

The changes in file format and the VCF convert tool have made the tool more usable for a broader range of applications, and by calculating the complete
two-tailed FET and implementing three different methods for calculating MDS, we have achieved a more accurate result. 
We found that we were able to gain a large speedup of our tools, giving a dramatic decrease in run time. 
This should make it possible to run several analyses with different parameters in a short amount of time. The tools were able
to reproduce the results found by previous studies, and produce good results on two additional data sets. The CSS was the most accurate method, achieving
good results over a broad range of data sets, while the FET had more noise and is weaker on data sets with little genomic divergence.

The increased speed and user friendliness of the tools make it feasible to run these analyses on a larger scale than has previously been done. These tools should be
able to meet the increased need for analysing large scale data sets in comparative genomics.
