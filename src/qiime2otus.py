#!/usr/bin/env python

import csv
import os
import re
import sys

c_CapMemoryUse		= 10,000 #lines of data that will be saved between sending data to the file
c_strUnclassified	= "unclassified"
c_strOTUID		= "#OTU ID"
c_cOutputOTULineageDelim = "|"
c_cOutputOTULineageDelimForRe = "\|"
c_strOutputDelimiter	="\t"

#Extensions
c_strQiime1Extension = ".qiime"
c_strQiime2Extension = ".qiime2"

#Qiime Format 1
def qiimeFormat1(strConsensusLineage):
  #Format consensus lineage
  #Remove root
  strConsensusLineage = re.sub(r'^$',c_strUnclassified,strConsensusLineage)
  #Remove root
  strConsensusLineage = re.sub(r'^',"",strConsensusLineage)
  #Change no end clade to unclassified
  strConsensusLineage = re.sub(r'(;)+$',c_cOutputOTULineageDelim+c_strUnclassified,strConsensusLineage)
  #Change otu delimiters
  strConsensusLineage = re.sub(r';',c_cOutputOTULineageDelim,strConsensusLineage)
  #Indicate internal unclassifieds
  strNewLineage = ""
  while(not strNewLineage == strConsensusLineage):
    strNewLineage = strConsensusLineage
    strConsensusLineage = re.sub(r'\|\|',c_cOutputOTULineageDelim+c_strUnclassified+c_cOutputOTULineageDelim,strConsensusLineage)
  #If the consensus lineage was only an inital root make sure it is called unclassified and not blank
  if strConsensusLineage == "":
    strConsensusLineage = c_strUnclassified
  return strConsensusLineage

#Qiime Format 2
def qiimeFormat2(strConsensusLineage):
  #Format consensus lineage
  #Remove root
  strConsensusLineage = re.sub(r'^Root$',c_strUnclassified,strConsensusLineage)
  #Remove root
  strConsensusLineage = re.sub(r'^Root;',"",strConsensusLineage)
  #Change no end clade to unclassified
  strConsensusLineage = re.sub(r'(;)+$',c_cOutputOTULineageDelim+c_strUnclassified,strConsensusLineage)
  #Change otu delimiters
  strConsensusLineage = re.sub(r';',c_cOutputOTULineageDelim,strConsensusLineage)
  #Indicate internal unclassifieds
  strNewLineage = ""
  while(not strNewLineage == strConsensusLineage):
    strNewLineage = strConsensusLineage
    strConsensusLineage = re.sub(r'\|\|',c_cOutputOTULineageDelim+c_strUnclassified+c_cOutputOTULineageDelim,strConsensusLineage)
  #If the consensus lineage was only an inital root make sure it is called unclassified and not blank
  if strConsensusLineage == "":
    strConsensusLineage = c_strUnclassified
  return strConsensusLineage

#Check for arg count
if len( sys.argv ) != 3:
	raise Exception( "Usage: qiime2otus.py <inputfile.txt> <outputfile.otu>" )
strInputQiime = sys.argv[1]
strOutputName = sys.argv[2]

print("Outputing data to the file:"+strOutputName)

#Temporarly holds lines of data to output
lstrOutputlines = []

#If output file exists, erase
if os.path.exists(strOutputName):
  os.remove(strOutputName)

#Get file extension
sFileExtension = "".join([".",filter(None,strInputQiime.split("."))[-1]])

#Tracks if errors occured
sErrors = 0

#Output file handle
fhndlOutput = open(strOutputName, "a")

#Read through tab delimited qiime output
for astrLine in csv.reader( open(strInputQiime), csv.excel_tab ):
  #Get otu and data
  strOTU, astrData, strConsensusLineage = astrLine[0], astrLine[1:-1], astrLine[-1]

  #Consolidate the data
  astrData = c_strOutputDelimiter.join(astrData)

  #If the OTU id is found save, if data is found save, ignore other header elements
  if strOTU[0] == "#":
    if strOTU == c_strOTUID:
      lstrOutputlines.append(c_strOutputDelimiter.join([strOTU,astrData]))
  #If the otu name has a period in it, take everything after the period
  else:
    i = strOTU.find( "." )
    if i >= 0:
      strOTU = strOTU[( i + 1 ):]
    #Format consensus lineage based on extension
    if sFileExtension == c_strQiime1Extension:
      strConsensusLineage = qiimeFormat1(strConsensusLineage)
    elif sFileExtension == c_strQiime2Extension:
      strConsensusLineage = qiimeFormat2(strConsensusLineage)
    else:
      print("Error unknown file extension:"+sFileExtension+". Was not translated.")
      sErrors = sErrors + 1
      break
    #Combine with qiime consensus lineages
    strOTU = c_cOutputOTULineageDelim.join([strConsensusLineage,strOTU])
    lstrOutputlines.append(c_strOutputDelimiter.join([strOTU,astrData]))

  #If the cache of data is a certain size then output to file
  if len(lstrOutputlines)>c_CapMemoryUse:
    fhndlOutput.write("\n".join(lstrOutputlines))
    lstrOutputlines = []

#Output to file anthing not already sent to file
fhndlOutput.write("\n".join(lstrOutputlines))
fhndlOutput.close()

if sErrors:
  print("Errors occured during Qiime2OTU.")
else:
  print("Successfully Ended Qiime2OTU.")
