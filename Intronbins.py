#!/usr/bin/env python3

#Import modules
import sys
import pandas
import csv

#Initialise variables / bins
#1 = 1-20
Bin1=0
#2 = 21-50
Bin2=0
#3 = 51-100
Bin3=0
#4 = 101-200
Bin4=0
#5 = 201-400
Bin5=0
#6 = >400
Bin6=0
count=0
cumFeatLength=0

#Import and parse the file
input = sys.argv[1]
count = len(open(input).readlines(  ))

lol = list(csv.reader(open(input, 'rt'), delimiter='\t'))
for i in range(0,count):
	#TE locations
	LTE=lol[i][3]
	RTE=lol[i][4]
	TElength=int(RTE)-int(LTE)
	#print(LTE)
	#print(RTE)	

	#TE name
	TE=lol[i][8]
	#Processing (due to not being tab delimited)
	TEsplit1 = TE.split("\"")
	TEsplit2 = TEsplit1[1].split(":")
	#print(TEsplit2[1])

	#Feature locations
	LFeat=lol[i][12]
	RFeat=lol[i][13]

	#Extracting intron value
	Intron=lol[i][17]
	Isplit1 = Intron.split(";")
	Isplit2 = Isplit1[0].split("=")
	Isplit3 = Isplit1[1].split("=")	
	Intront = Isplit2[1].split("_added-")
	IntronID = str(Intront[0]) + "-" +str(Intront[1])

	#Extracting parent value, just incase?
	Parent = Isplit3[1]

	#Feature Lengths
	FeatLength=abs(int(LFeat)-int(RFeat))
	cumFeatLength+=FeatLength

	#Sanity Checks
	#print(FeatLength)
	#print(IntronID)
	#print(Parent)
	#print(LFeat)
	#print(RFeat)
	#print(str(IntronID) +  " is " +  str(FeatLength) + " base pairs long and is in " + Parent)	
	#print(str(TEsplit2[1]) +  " fragment is " +  str(TElength) + " base pairs long")


	#Distance from coding regions / work out direction
	
	#Check which side is closest to what thing
	Range1=abs(int(LFeat)-int(LTE))
	Range2=abs(int(LFeat)-int(RTE))
	Range3=abs(int(RFeat)-int(LTE))
	Range4=abs(int(RFeat)-int(RTE))
	#print(Range1)
	#print(Range2)
	#print(Range3)
	#print(Range4)
	
	if(Range1<Range3):
		#LFeat & LTE
		distance=Range1
	else:
		#RFeat & RTE
		distance=Range4

	#distance = bp to closest coding region

	#Sanity Checks 2	
	#print(distance)
	#print(LFeat)
	#print(LTE)
	#print(RFeat)
	#print(RTE)

	#Range2/4 can be used to ascertain orientation of element
	#if(Range1<Range2):
		#print("Normal Orientation")


	#Add to bins
	#1 = 1-20, 2 = 21-50, 3 = 51-100, 4 = 101-200, 5 = 201-400, 6 = >400
	if(distance<21):
		Bin1+=1
	elif(distance<51):
		Bin2+=1
	elif(distance<101):
		Bin3+=1
	elif(distance<201):
		Bin4+=1
	elif(distance<401):
		Bin5+=1
	else:
		Bin6+=1

#Present results
print("Bin 1, 1-20bp: "+ str(Bin1))
print("Bin 2, 21-50bp: "+ str(Bin2))
print("Bin 3, 51-100bp: "+ str(Bin3))
print("Bin 4, 101-200bp: "+ str(Bin4))
print("Bin 5, 201-400bp: "+ str(Bin5))
print("Bin 6, >400bp: "+ str(Bin6))
	
StdBin1 = Bin1/20
StdBin2 = Bin2/30
StdBin3 = Bin3/50
StdBin4 = Bin4/100
StdBin5 = Bin5/200
#Defo wrong this last one
StdBin6 = Bin6/9271
#print(cumFeatLength)
AvgFeatLength=cumFeatLength/count
#print(AvgFeatLength)

print("Standardised")
print("StdBin 1, 1-20bp: "+ str(StdBin1))
print("StdBin 2, 21-50bp: "+ str(StdBin2))
print("StdBin 3, 51-100bp: "+ str(StdBin3))
print("StdBin 4, 101-200bp: "+ str(StdBin4))
print("StdBin 5, 201-400bp: "+ str(StdBin5))
print("StdBin 6, >400bp: "+ str(StdBin6))








