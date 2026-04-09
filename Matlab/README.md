The main folder here contains MATLAB scripts for SNV, LSNV, PSNV, MSC, LMSC, and PMSC, along with a readme and demo file. 

The 'Other Scripts' folder contains some alternative versions of PMSC and PSNV, including using C++ code. These will give identical results (within machine tolerance) to the main pmsc.m and psvn.m scripts, but are less computationally efficient in MATLAB. 
psnv.m uses the built-in MATLAB functions 'movmean' and 'movstd' for rapid computation. pmsc.m performs moving PMSC on a sliding band using vectorization and achieves fast results by avoiding the need for loops. 
