
----------------------------------------------------------------------

LOAD A DATA SET (requires the DataSet Object (DSO), available at https://www.mathworks.com/matlabcentral/fileexchange/39336-dataset-object) 

>> corn = load('corn.mat');

----------------------------------------------------------------------

YOU SHOULD SEE THE FOLLOWING.

>> corn

corn = 

  struct with fields:

    information: [18×60 char]
         m5spec: [80×700 dataset]
        mp5spec: [80×700 dataset]
        mp6spec: [80×700 dataset]
       propvals: [80×4 dataset]
          m5nbs: [3×700 dataset]
         mp5nbs: [4×700 dataset]
         mp6nbs: [4×700 dataset]

----------------------------------------------------------------------

LOAD THE SPECTRA (MP5) AND WAVELENGTHS

>> X = corn.mp5spec.data;
>> wave = corn.mp5spec.axisscale{2,1};

----------------------------------------------------------------------

PLUG THIS DATA INTO THE DEMO CODE.  THE DEMO CODE RUNS ALL OF THE 
MSC AND SNV VARIANTS.  THE DEMO CODE GIVES EXAMPLES ON THE EXPECTED 
INPUT FOR EACH FUNCTION.  A FIGURE WITH THREE SUBPLOTS WILL BE 
GENERATED.  
[1] THE UPPER LEFT: ORIGINAL SPECTRA
[2] THE LOWER LEFT: TRANSFORMED SPECTRA ASSOCIATED WITH 
    HELLAND-BASED MSC AND ITS VARIANTS (HELLAND-BASED MSC AND 
	CLASSICAL MSC ARE EQUIVALENT)
[3]	THE LOWER RIGHT: TRANSFORMED SPECTRA ASSOCIATED WITH 
    SNV AND VARIANTS.

>> demo(X,wave)

NOTE: THE FIGURE GENERATED WILL LIKELY HAVE TO BE ENLARGED TO CLEARLY
THE DIFFERENCES ACROSS SUBPLOTS.
