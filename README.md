# WignerPR: Matlab Software for Phase Retrieval using the Wigner Deconvolution Method

This repository contains Matlab code for solving the (ptychographic) 
phase retrieval problem using the Wigner distribution deconvolution 
method. Details of the method, theoretical guarantees and 
representative numerical results can be found in the following 
manuscript      

'Inverting Spectrogram Measurements via Aliased Wigner Distribution
Deconvolution and Angular Synchronization' 
Michael Perlmutter, Sami Merhi, Aditya Viswanthan and Mark Iwen
https://arxiv.org/abs/1907.10773


This software was developed at the [Department of 
Mathematics][msumath], [Michigan State University][msu] 
and the [Department of Mathematics and Statistics][uofmdmath], 
[University of Michigan - Dearborn][uofmd] and is released under 
the MIT license.

The software was developed and tested using Matlab 
R2019a. 


## Directory Structure and Contents

The software package is organized under the following 
directory structure:

 - figures/    
   This folder (and its subfolders) contains Matlab scripts to 
   regenerate plots from the manuscript.

 - supplementary_numerics/    
   This folder contains Matlab scripts and additional plots 
   and documents which summarize further numerical results of 
   potential interest to the reader. (Consult the ReadMe.pdf 
   in the folder for further details.)

 - src/    
   This folder contains (self-contained) simple Matlab 
   implementations of the methods proposed in the manuscript.


## Instructions

Extract the contents of the zip file and execute, in 
Matlab/Octave, scripts from the folders. 


## Contact

Bug reports, comments and suggestions are welcome 
at the [BlockPR Bitbucket repository][bitbucket]. 


[msu]: http://www.msu.edu/
[msumath]: http://math.msu.edu/
[uofmd]: https://umdearborn.edu/
[uofmdmath]: https://umdearborn.edu/casl/departments/mathematics-and-statistics
[bitbucket]: https://bitbucket.org/charms/blockpr/
