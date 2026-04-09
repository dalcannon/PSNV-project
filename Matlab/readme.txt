
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

>> X = mp5spec.data;
>> wave = mp5spec.axisscale{2,1};

----------------------------------------------------------------------

PLUG THIS DATA INTO THE DEMO CODE.  THE DEMO CODE RUNS ALL OF THE 
MSC AND SNV VARIANTS.  THE DEMO CODE GIVES EXAMPLES ON THE EXPECTED 
INPUT FOR EACH FUNCTION.  A FIGURE WITH FOUR SUBPLOTS WILL BE 
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
