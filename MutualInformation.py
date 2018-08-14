#This code does: 
#1) discretizes continuous data by binning using linspace
#2) adds a column of 0's and 1's to uninduced and induced, respectively; this allows calculation of the entropy of the input, "input"
#3) combines the data for uninduced, x, and induced, y, into "mergeddata"
#4) calculates entropy using the formula entropyx = -sum(p(x)*log2 p(x))
#5) calculates mutual information using the formula entropy mergeddata + entropy input - [entropy (x)+entropy(y)]
#6) Note: this processes files with the name "Date_StrainName_LineNumber_PassageNumber_Condition_x_x_x"; 

import numpy as np
import matplotlib
#import matplotlib.pyplot as plt #can't import this when running on terminal (machine without display)
import csv
import glob
import os
import datetime #to save outfile with today's date

matplotlib.use('Agg') #calling this to fix an error with opening on machines without display

import matplotlib.pyplot as plt #calling after previous line to fix an error with opening on machines without display
import matplotlib.mlab as mlab
from matplotlib.patches import Rectangle

#Use FindFewestRowsLowestandHighestGFP.R to find lower and upper limits for GFP values
lowerlimit=0
upperlimit=5.5
#Number of bins should be at least ~0.05% of the number of data points. Example: at least 1000 bins for 20,000 cells
numberofbins="1000"
path = "" #insert path name here

now = datetime.datetime.now()
todaysdate = now.strftime("%Y%m%d") #get's today's date in the format yyyymmdd

outfilename = "MutualInformation_"+todaysdate+"_Range" + str(lowerlimit) + "-" + str(upperlimit) + "_" + str(numberofbins) + "bins.csv"
outfile = open(outfilename, "w")  # opens and with "w", writes, a new file in the same directory as path with the given name
# write only takes one parameter so adding all the different column names together 
# (don't need to here but just doing it for consistency with variables later)
# separating with a "," so the file can be csv
# printing a new line at the end so the values can be propagated in a new row
outfile.write('Sample1,' + 'Sample2,' + 'MI_GFP,' + 'Date,' + 'Strain,' + 'Line,' + 'Passage,' + 'Time,\n') #add MI_BFP if necessary

#Calculates entropy of the merged data 
def calculate_entropy (data, cells, name, lowerlimit, upperlimit, numberofbins):

	bins = np.linspace(lowerlimit, upperlimit, numberofbins) #bins values using the limits into the given number of bins
	digitized = np.digitize(data, bins) #Return the indices of the bins to which each value in input array belongs
	counts = np.bincount(digitized)[1:].astype(np.float32) #Count number of occurrences of each value in array of non-negative ints.
	normalized = counts / cells #calculates p(x) 

	sum = 0
	for i in normalized:
		if i != 0:
			sum += i * np.log2(i)
		elif i == 0: 
			sum += i

	if name != 'none':
		name_split = name.split("_")
		date=name_split[0]
		strain=name_split[1]
		line=name_split[2]
		passage=name_split[3]
		#condition=name_split[4]
		time=name_split[5]
		time=time+"min"
		reporter=name_split[12]
		figuretitle = date+" "+strain+" "+line+" "+passage+" "+reporter+" response \n+&- Induction after "+time+" ("+numberofbins+" bins)"

		plt.hist(data, bins=bins, color='red', alpha=0.5, edgecolor='red')
		plt.title(figuretitle)
		plt.ylabel("Count")
		plt.xlabel("Fluorescent Reporter Expression (A.U.)")

	return -sum

#Calculates entropy of the input
def calculate_entropy_input (data, cells):
	bins = np.linspace(0, 1, 2) 
	digitized = np.digitize(data, bins)
	counts = np.bincount(digitized)[1:].astype(np.float32)

	normalized = counts / cells

	sum = 0
	for i in normalized:
		if i != 0:
			sum += i * np.log2(i)
		elif i == 0: 
			sum += i

	return -sum

#Calculates joint entropy
def calculate_entropy_xy (data_x, data_y, cells, xname, yname, lowerlimit, upperlimit, numberofbins):

	bins_x = np.linspace(lowerlimit, upperlimit, numberofbins) 
	digitized_x = np.digitize(data_x, bins_x) 
	counts_x = np.bincount(digitized_x)[1:].astype(np.float32) 

	bins_y = np.linspace(lowerlimit, upperlimit, numberofbins) 
	digitized_y = np.digitize(data_y, bins_y)
	counts_y = np.bincount(digitized_y)[1:].astype(np.float32)

	normalized_x = counts_x / cells
	normalized_y = counts_y / cells

	sum = 0
	sum_x = 0
	for i in normalized_x:
		if i != 0:
			sum_x += i * np.log2(i)
		elif i == 0:
			sum_x += i

	sum_y = 0
	for i in normalized_y:
		if i != 0:
			sum_y += i * np.log2(i)
		elif i == 0:
			sum_y += i

	sum = sum_x+sum_y

	name = yname
	name_split = name.split("_")
	date = name_split[0]
	strain = name_split[1]
	line = name_split[2]
	passage = name_split[3]
	#condition = name_split[4]
	time = name_split[5]
	time = time + "min"
	reporter = name_split[6]

	#plt.figure()
	plt.hist(data_x, bins=bins_x, color='grey', alpha=0.8, edgecolor='grey')
	#plt.title(figuretitlex)
	#plt.savefig(xname+".png")
	#plt.close()

	#plt.figure()
	plt.hist(data_y, bins=bins_y, color='black', alpha=0.7)
	#plt.title(figuretitle)
	handles = [Rectangle((0, 0), 1, 1, color=c, ec="k") for c in ['grey', 'black', 'red']]
	labels = ["Uninduced", "Induced", "Combined"]
	plt.legend(handles, labels, loc=0)
	now = datetime.datetime.now()
	todaysdate = now.strftime("%Y%m%d")
	figurename = "Range"+str(lowerlimit)+"-"+str(upperlimit)+"_"+numberofbins + "bins_"+date + "_" + strain + "_" + line + "_"+ passage + "_" + reporter +"_"+ time+"_"+todaysdate+".png"
	plt.savefig(figurename)
	plt.close()

	return -sum

def open_file (path):
    file1 = open(path)
    list1 = csv.reader(file1, delimiter=',')
    matrix = np.array([row for row in list1])
    return matrix

#copied this function from: https://mail.python.org/pipermail/tutor/2004-April/029019.html
#Mac OS X inserts .DS_Store files and this function ignores any files that begin with a "."
def mylistdir(directory):
    """A specialized version of os.listdir() that ignores files that
    start with a leading period."""
    filelist = os.listdir(directory)
    return [x for x in filelist
            if not (x.startswith('.'))]

for filename in mylistdir(path):
	print("Filename: ", filename)
	strainName = None
	splitResult = None
	date = None
	line = None
	passage = None
	strain = None
	condition = None
	time = None
	rest = None
  	path1name = None
	path1 = None
	path2name = None
	path2 = None

	if filename.endswith(".csv"):
		strainName = filename
		splitResult = filename.split("_")  # separates the filename whereever there's an underscore
		date = splitResult[0]
		strain = splitResult[1]
		line = splitResult[2]
		passage = splitResult[3]
		condition = splitResult[4]
		time = splitResult[5]
		rest = ' '.join(["_",splitResult[6],"_",splitResult[7],"_",splitResult[8]])

		if condition == 'Uninduced':
			path1name = ' '.join([date, "_", strain, "_", line, "_",passage,"_Uninduced", "_", time, rest])
			path1name = path1name.replace(" ", "")  # remove spaces
			path1 = path1name
			path2name = ' '.join([date, "_", strain, "_", line, "_",passage,"_Induced", "_", time, rest])  # add the string together
			path2name = path2name.replace(" ", "")  # remove spaces
			path2 = path2name

		else:
			continue

		print("Path1: ", path1)
		data1 = open_file(path1)
		print("Path2: ", path2)
		data2 = open_file(path2)

		# Mutual Information (GFP)
    #Stores GFP values for uninduced cells from column 7 into a new array and adds a column of 0s indicating induction status=uninduced
		data1_GFP = data1[1:, 6]  # gets GFP "B1-A" column, count starts at 0 so it's column 6
		data1_GFP_noblanks = list(filter(None, data1_GFP))  # remove blanks
		data1_GFP_noblanks_array = np.asarray(data1_GFP_noblanks)  # convert list to array
		data1_GFP_noblanks_array = data1_GFP_noblanks_array.astype(np.float32)
		zeros_GFP = np.zeros((len(data1_GFP_noblanks_array), 1))
		data1_GFP_noblanks_array_reshape = np.reshape(data1_GFP_noblanks_array, (len(data1_GFP_noblanks_array), 1))  # changes the shape of the array
		data1_GFP_noblanks_array_condition = np.append(data1_GFP_noblanks_array_reshape, zeros_GFP, axis=-1)
    
    #Stores GFP values for induced cells from column 7 into a new array and adds a column of 1s indicating induction status=induced
		data2_GFP = data2[1:, 6]
		data2_GFP_noblanks = list(filter(None, data2_GFP))
		data2_GFP_noblanks_array = np.asarray(data2_GFP_noblanks)
		data2_GFP_noblanks_array = data2_GFP_noblanks_array.astype(np.float32)
		ones_GFP = np.ones((len(data2_GFP_noblanks_array), 1))
		data2_GFP_noblanks_array_reshape = np.reshape(data2_GFP_noblanks_array, (len(data2_GFP_noblanks_array), 1))
		data2_GFP_noblanks_array_condition = np.append(data2_GFP_noblanks_array_reshape, ones_GFP, axis=-1)
    
    #Merging uninduced and induced
		mergeddata_GFP = np.concatenate((data1_GFP_noblanks_array_condition, data2_GFP_noblanks_array_condition))
		mergeddata_GFP = np.array(mergeddata_GFP)
		totalcells_GFP = len(mergeddata_GFP) #gets length of merged data array to find total number of cells
    
    #Gets file names without full path
		data1name = os.path.basename(path1)
		data1name = data1name.replace(".csv", "")
		data2name = os.path.basename(path2)
		data2name = data2name.replace(".csv", "")
		mergedname = data1name + "_" + data2name
    
    #Calculate entropy
		entropy1_GFP = calculate_entropy(mergeddata_GFP[:, 0], totalcells_GFP, (mergedname + "_GFP"), lowerlimit, upperlimit, numberofbins)
		entropy2_GFP = calculate_entropy_input(mergeddata_GFP[:, 1], totalcells_GFP)
		joint_entropy_GFP = calculate_entropy_xy(data1_GFP_noblanks_array, data2_GFP_noblanks_array, totalcells_GFP, (data1name + "_GFP"), (data2name + "_GFP"), lowerlimit, upperlimit, numberofbins)
    
    #Calculate MI
		mutual_information_GFP = entropy1_GFP + entropy2_GFP - joint_entropy_GFP
		mutual_information_GFP = str(mutual_information_GFP)  # convert number to string so it can be stored (issue with storing mixed number and string array)

#		#Mutual Information (BFP)
    #Stores BFP values for uninduced cells from column 8 into a new array and adds a column of 0s indicating induction status=uninduced
#		data1_BFP = data1[1:, 7] 
#		data1_BFP_noblanks = list(filter(None, data1_BFP))  # remove blanks
#		data1_BFP_noblanks_array = np.asarray(data1_BFP_noblanks)  # convert list to array
#		data1_BFP_noblanks_array = data1_BFP_noblanks_array.astype(np.float32)
#		zeros_BFP = np.zeros((len(data1_BFP_noblanks_array), 1))
#		data1_BFP_noblanks_array_reshape = np.reshape(data1_BFP_noblanks_array, (len(data1_BFP_noblanks_array), 1))  # changes the shape of the array
#		data1_BFP_noblanks_array_condition = np.append(data1_BFP_noblanks_array_reshape, zeros_BFP, axis=-1)
#
#   #Stores BFP values for induced cells from column 8 into a new array and adds a column of 1s indicating induction status=induced
#		data2_BFP = data2[1:, 7]
#		data2_BFP_noblanks = list(filter(None, data2_BFP))
#		data2_BFP_noblanks_array = np.asarray(data2_BFP_noblanks)
#		data2_BFP_noblanks_array = data2_BFP_noblanks_array.astype(np.float32)
#		ones_BFP = np.ones((len(data2_BFP_noblanks_array), 1))
#		data2_BFP_noblanks_array_reshape = np.reshape(data2_BFP_noblanks_array, (len(data2_BFP_noblanks_array), 1))
#		data2_BFP_noblanks_array_condition = np.append(data2_BFP_noblanks_array_reshape, ones_BFP, axis=-1)
#   
#   #Merging uninduced and induced
#		mergeddata_BFP = np.concatenate((data1_BFP_noblanks_array_condition, data2_BFP_noblanks_array_condition))
#		mergeddata_BFP = np.array(mergeddata_BFP)
#		totalcells_BFP = len(mergeddata_BFP)
# 
#   #Calculate Entropy  
#		entropy1_BFP = calculate_entropy(mergeddata_BFP[:, 0], totalcells_BFP, (mergedname + "_BFP"), lowerlimit, upperlimit, numberofbins)
#		entropy2_BFP = calculate_entropy_input(mergeddata_BFP[:, 1], totalcells_BFP)
#		joint_entropy_BFP = calculate_entropy_xy(data1_BFP_noblanks_array, data2_BFP_noblanks_array, totalcells_BFP, (data1name + "_BFP"), (data2name + "_BFP"), lowerlimit, upperlimit, numberofbins)
#   #Calculate Mutual Information
#		mutual_information_BFP = entropy1_BFP + entropy2_BFP - joint_entropy_BFP
#		mutual_information_BFP = str(mutual_information_BFP)  # convert number to string so it can be stored (issue with storing mixed number and string array)

		# for testing:
		# print ('\tPath1:', os.path.basename(path1))
		# print ('\tPath2:', path2)
		# print ('\t\tmutual_information_GFP:', mutual_information_GFP)
		# print ('\t\tmutual_information_BFP:', mutual_information_BFP)
		# print ()
		#conditionname = condition.replace("No", "")
		
		#add mutual_information_BFP to the output if necessary
    		outfile.write(os.path.basename(path1) + ',' + path2 + ',' + mutual_information_GFP +  ',' + date + ',' + strain + ',' + line + ',' + passage +',' + time + ',\n') 
		plt.close('all')

outfile.close()
