Add some real datasets, and some Matlab/R scripts with the results from their analysis

A brief summary of the datasets relevant to the paper: 

Marzipan dataset: 32 samples of marzipan. Pretty impressive prediction performance for regional methods, but low sample size

Tablet dataset: 310 tablets. Most of the chemical info from the API is at a small region of the spectrum, and so regional scatter correction methods will have an advantage in these cases. 
Prediction performance is pretty good for regional scatter correction methods. 
Most difficult part was that it's got some Experimental Design aspects, and it wasn't super obvious which segments (Pilot-scale, full-scale, lab-scale) should go in training set vs. test set. 
I'll probably need to slightly redo my results in this section, but don't expect RMSEP's to change very much. 

Tecator Data: 215* samples, but only 100 wavelengths. Works pretty well for predicting Moisture & Fat, but not Protein. 
*10% of the 215 samples are exact copies of other samples. 
