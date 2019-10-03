Change in Mutual Information Calculator
===
Author: Purnima S. Kompella purnima.kompella@mail.utoronto.ca

To install in unix based systems:

```
git clone https://github.com/purnimakompella/MutualInformation.git
```

Step 1: Convert fcs files to csv
1) In FCStolog10CSV_Fortessa.R, insert path name of folder containing fcs files in line 11: parent.dir <- ""
Note: This program gets values from columns labeled "FSC-A", "V_450/50-A", "B_530/30-A". 

Step 2: Obtain minimum cell count, and range of fluorescent reporter values
Note: The column with reporter values in the csv files should be labeled 'B_530_30.A'
1) In FindFewestRowsLowestandHighestGFP.R, insert path name of folder containing CSV files in line 9: parent.dir <- ""
2) Run the program. This will print minimum number of rows in all the files, and the minimum and maximum fluorescent reporter values.

Step 3: Calculate change in mutual information
In mutualinformation.py:
1) Replace the 0 with minimum reporter value in line 24: lowerlimit=0
2) Replace the 5.5 with maximum reporter value in line 25: upperlimit=5.5
3) Replace 1000 with atleast 10% of minimum number of rows in line 26: numberofbins="1000" [Example: at least 1000 bins for 20,000 cells]
4) Run mutualinformation.py with path name of file containing fcs files as an argument [Example: python mutualinformation.py pathname]

This code does the following: 
1) discretizes continuous data by binning using linspace
2) adds a column of 0's and 1's to uninduced and induced, respectively; this allows calculation of the entropy of the input, "input"
3) combines the data for uninduced, x, and induced, y, into "mergeddata"
4) calculates entropy using the formula entropyx = -sum(p(x)*log2 p(x))
5) calculates mutual information using the formula entropy mergeddata + entropy input - [entropy (x)+entropy(y)]
6) Note: this processes files with the name "Date_StrainName_LineNumber_PassageNumber_Condition_x_x_x"; 


