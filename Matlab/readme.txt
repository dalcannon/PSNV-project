
----------------------------------------------------------------------

LOAD A DATA SET

>> load corn.mat;

----------------------------------------------------------------------

YOU SHOULD SEE THE FOLLOWING.

>> corn

corn = 

  struct with fields:

           Xm5: [80×700 double]
          Xmp5: [80×700 double]
          Xmp6: [80×700 double]
       Xm5_nbs: [3×700 double]
      Xmp5_nbs: [4×700 double]
      Xmp6_nbs: [4×700 double]
    attributes: {'moisture'  'oil'  'protein'  'starch'}
             y: [80×4 double]
          wave: [1100 1102 1104 1106 1108 1110 … ] (1×700 double)
          info: [7×59 char]

----------------------------------------------------------------------

LOAD THE SPECTRA (MP5) AND WAVELENGTHS

>> X = corn.Xm5;
>> wave = corn.wave;

----------------------------------------------------------------------

PLUG THIS DATA INTO THE DEMO CODE.  THE DEMO CODE RUNS ALL OF THE 
MSC AND SNV VARIANTS.  THE DEMO CODE GIVES EXAMPLES ON THE EXPECTED 
INPUT FOR EACH FUNCTION.  A FIGURE WITH FOUR SUBPLOTS WILL BE 
GENERATED.  
[1] THE UPPER LEFT: ORIGINAL SPECTRA
[2] THE UPPER RIGHT: TRANSFORMED SPECTRA ASSOCIATED WITH 
    CLASSICAL MSC AND ITS VARIANTS
[3] THE LOWER LEFT: TRANSFORMED SPECTRA ASSOCIATED WITH 
    HELLAND-BASED MSC AND ITS VARIANTS (HELLAND-BASED MSC AND 
	CLASSICAL MSC ARE EQUIVALENT)
[4]	THE LOWER RIGHT: TRANSFORMED SPECTRA ASSOCIATED WITH 
    SNV AND VARIANTS.

>> demo(X,wave)

NOTE: THE FIGURE GENERATED WILL LIKELY HAVE TO BE ENLARGED TO CLEARLY
THE DIFFERENCES ACROSS SUBPLOTS.

----------------------------------------------------------------------

MATLAB SCRIPTS WITH THE SUFFIX _v2.m CORRESPOND TO CALLS TO EXTERNAL 
MEX BINARY FILES ENDING WITH THE SUFFIX _v2mex.mexw64.  THIS BINARY 
WAS COMPILED ON A WINDOWS 64-BIT MACHINE.  IF YOU HAVE LINUX OR MAC 
OS, THEN THE CORRESPONDING C++ FILE (WITH SUFFX v2mex.cpp) IS 
PROVIDED SO THAT YOU CAN COMPILE THE FILE ON YOUR OPERATING SYSTEM, 
I.E.,

>> mex pmsc_ccc_v2mex.cpp

