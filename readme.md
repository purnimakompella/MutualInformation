Mutual Information Project
===
Author: Purnima S. Kompella

This code does: 
1) discretizes continuous data by binning using linspace
2) adds a column of 0's and 1's to uninduced and induced, respectively; this allows calculation of the entropy of the input, "input"
3) combines the data for uninduced, x, and induced, y, into "mergeddata"
4) calculates entropy using the formula entropyx = -sum(p(x)*log2 p(x))
5) calculates mutual information using the formula entropy mergeddata + entropy input - [entropy (x)+entropy(y)]
6) Note: this processes files with the name "Date_StrainName_LineNumber_PassageNumber_Condition_x_x_x"; 

change