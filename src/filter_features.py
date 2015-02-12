#!/usr/bin/env python

import csv
import numpy
import sys

if len( sys.argv ) != 4:
	raise Exception( "Usage: filter_features.py <metadatum> <mean> <max> < <otus.txt>" )
strMetadatum, strMean, strMax = sys.argv[1:]
dMean, dMax = (float(s) for s in (strMean, strMax))

csvw = csv.writer( sys.stdout, csv.excel_tab )
fData = False
for astrLine in csv.reader( sys.stdin, csv.excel_tab ):
	strID, astrData = astrLine[0], astrLine[1:]
	if fData:
		adData = [float(s) for s in astrData]
		if ( max( adData ) < dMax ) or ( numpy.mean( adData ) < dMean ):
			continue
	elif strID == strMetadatum:
		fData = True
	csvw.writerow( [strID] + astrData )
