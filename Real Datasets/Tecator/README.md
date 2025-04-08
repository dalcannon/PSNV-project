Tecator meat dataset (https://lib.stat.cmu.edu/datasets/tecator)
tecator_X.csv and tecator_y.csv contain the original file with 215 spectra. 
tecax gives the wavelengths

However, the original dataset has duplicate samples which are exact copies of one another (for both X & y). These pairs of samples are as follows: 
Duplicates: 
        1  2  3  4  5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22
Dupe1: 28 12 15 16 61  36  37  38  40  58  75  26  39  64 145  86  89  13  17  59 139 186
Dupe2: 29 48 51 54 62 110 116 117 123 135 136 138 142 150 176 180 181 188 190 192 204 215

I also have attached tecXrm and tecyrm, which are versions of these datasets with one member of each duplicate pair removed (189 samples). 
