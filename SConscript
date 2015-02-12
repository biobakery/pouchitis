import csv
import os
import re
import sfle
import sys

Import( "*" )
pE = DefaultEnvironment( )

c_dMean							= 0.005
c_dMax							= 0.05
c_iMinCount						= 3
c_iMinSamples					= 3
c_strMetadatum					= "ileo close"
c_strKeys						= "Location|Outcome"
c_strREFeatures					= r'^((\d+)|(Unclassified.+)|(k__.+))$'
c_strPre						= "_pre"

c_fileInputOTUsTXT				= sfle.d( pE, fileDirInput, "otu_table.txt.gz" )
c_fileInputGenesPCL				= sfle.d( pE, fileDirInput, "U100-pouchitus_82012_trimmed-transposed-checked.pcl" )
c_fileInputGeneMetadataTXT		= sfle.d( pE, fileDirInput, "gene expression metadata.txt" )
c_fileInputMetadataTXT			= sfle.d( pE, fileDirInput, "Sample List Broad_pheno7_6_2012.txt" )
c_fileInputGenesTXT				= sfle.d( pE, fileDirInput, "representative_genes.txt" ) 
#------------------------------------------------------------------------------ 
c_strRiskMetadatum				= "initialSampleID"
c_fileInputRiskOTUsTXT			= sfle.d( pE, fileDirInput, "otu_table.s2c50p1000.txt" )
c_fileInputRiskMapTXT			= sfle.d( pE, fileDirInput, "MAP_risk.txt" )
c_fileInputRiskMetadataTXT		= sfle.d( pE, fileDirInput, "serology_and_expression_8112012_subset.txt" )

c_fileGenesTXT					= sfle.d( pE, fileDirTmp, sfle.rebase( c_fileInputGenesTXT ) )
c_fileMetadataPCL				= sfle.d( pE, fileDirTmp, sfle.rebase( c_fileInputMetadataTXT, True, sfle.c_strSufPCL ) )
c_fileMetadataTXT				= sfle.d( pE, fileDirTmp, "metadata.txt" )
c_fileMetadataTSV				= sfle.d( pE, fileDirTmp, sfle.rebase( c_fileMetadataTXT, True, sfle.c_strSufTSV ) )
c_fileOTUsPCL					= sfle.d( pE, fileDirTmp, sfle.rebase( c_fileInputOTUsTXT, sfle.c_strSufGZ, "" ) )
c_fileMAssocsTXT				= sfle.d( pE, fileDirOutput, "massocs.txt" )
c_dirMaAsLin					= Dir( sfle.d( fileDirTmp, "maaslin" ) )

c_fileProgOTUs2PCL				= sfle.d( pE, fileDirSrc, "otus2pcl.py" )
c_fileProgFindIDs				= sfle.d( pE, fileDirSrc, "find_ids.py" )
c_fileProgFilterOTUs			= sfle.d( pE, fileDirSrc, "filter_otus.py" )
c_fileProgFilterFeatures		= sfle.d( pE, fileDirSrc, "filter_features.py" )
c_fileProgFilterSamples			= sfle.d( pE, fileDirSrc, "filter_samples.py" )
c_fileProgFileR					= sfle.d( pE, fileDirSrc, "maaslin_pouchitis.R" )
c_fileProgMultivariateMetadata	= sfle.d( pE, fileDirSrc, "multivariate_metadata.R" )
c_fileProgMAssocs2Table			= sfle.d( pE, fileDirSrc, "massocs2table.py" )
c_fileProgPlotLevels			= sfle.d( pE, fileDirSrc, "plot_levels.py" )
c_fileProgMergeMetadata			= sfle.find( c_fileDirInput, "merge_metadata.py" )
#------------------------------------------------------------------------------ 
c_fileProgRiskFileR				= sfle.d( pE, fileDirSrc, "maaslin_risk.R" )
c_fileProgQiime2PCL				= sfle.find( c_fileDirInput, "qiime2pcl.py" )

def sed( strFrom, strTo, fUnwrap = False ):
	def sed_helper( target, source, env, strFrom = strFrom, strTo = strTo ):
		strT, astrSs = sfle.ts( target, source )
		with open( strT, "w" ) as ostm:
			if fUnwrap:
				for astrLine in csv.reader( open( astrSs[0] ), csv.excel_tab ):
					for strToken in astrLine:
						strToken = re.sub( strFrom, strTo, strToken )
						ostm.write( "%s\n" % strToken )
			else:
				for strLine in open( astrSs[0] ):
					ostm.write( re.sub( strFrom, strTo, strLine ) )
		return False
	return sed_helper

Command( c_fileGenesTXT, c_fileInputGenesTXT, sed( "Repesentative_B", "sample", True ) )

#sfle.spipe( pE, c_fileInputGenesTXT, r"sed -E 's/\t/\n/g' | sed 's/^Representative_B/sample/'", c_fileGenesTXT )
apProcs = [
	sfle.CProcessor( "fix",	None,	sfle.CCommand( r"sed 's/-MBIOME-/.MBIOME./g'", None, True, True ),
		sfle.CTarget( None, fileDirTmp ) ),
	sfle.CProcessor( "trn",	None,	sfle.CCommand( c_fileProgTranspose ) ),
	sfle.CProcessor( "flt",	None,	sfle.CCommand( c_fileProgGrepRows,
		[[c_fileGenesTXT]] ),
		sfle.CTarget( None, fileDirOutput ) ),
]
fileGenesPCL = sfle.CProcessor.pipeline( pE, apProcs, c_fileInputGenesPCL )[0]

sfle.pipe( pE, c_fileInputOTUsTXT, c_fileProgFilterOTUs, c_fileOTUsPCL, [c_iMinCount, c_iMinSamples] )

sfle.pipe( pE, c_fileInputGeneMetadataTXT, c_fileProgTranspose, c_fileMetadataPCL )
sfle.cmd( pE, c_fileProgMergeTables, c_fileMetadataTXT, ["-t", [c_fileMetadataPCL], [fileGenesPCL]] )
apProcs = [
	sfle.CProcessor( "ids",	None,	sfle.CCommand( c_fileProgFindIDs ) ),
	sfle.CProcessor( "trn",	None,	sfle.CCommand( c_fileProgTranspose ) ),
]
fileGeneMetadataTSV = sfle.CProcessor.pipeline( pE, apProcs, c_fileMetadataTXT )[0]

#------------------------------------------------------------------------------ 

fileMetadataTSV = str(fileGeneMetadataTSV).replace( sfle.c_strSufTXT, sfle.c_strSufTSV )
sfle.cmd( pE, c_fileProgMergeTables, fileMetadataTSV, [[c_fileInputMetadataTXT], [fileGeneMetadataTSV]] )
fileDataPrePCL = sfle.CProcessor.pipeline( pE,
	sfle.CProcessor( "mtd",	None,	sfle.CCommand( c_fileProgMergeMetadata, [[fileMetadataTSV], "-m 0"] ) ),
	c_fileOTUsPCL )[0]
fileDataPCL = sfle.CProcessor.pipeline( pE,
	sfle.CProcessor( "flt",	None,	sfle.CCommand( c_fileProgFilterFeatures, [c_strMetadatum, c_dMean, c_dMax] ),
		sfle.CTarget( None, fileDirOutput ) ),
	fileDataPrePCL )[0]

sfle.pipe( pE, fileDataPrePCL, c_fileProgFilterFeatures,
	sfle.d( fileDirOutput, sfle.rebase( fileDataPrePCL, True, sfle.c_strSufPCL ) ),
	[c_strMetadatum, c_dMean / 10, c_dMax / 10] )

#------------------------------------------------------------------------------ 

def maaslin( dirMaAsLin, fileInputPCL, fileDataPCL, strFeatures, fileProgFileR, strMetadatum ):
	
	fileOTUsPreTSV = sfle.d( dirMaAsLin, sfle.rebase( fileInputPCL, True, c_strPre + sfle.c_strSufTSV ) )
	sfle.pipe( pE, fileDataPCL, c_fileProgTranspose, fileOTUsPreTSV )
	fileOTUsTSV = sfle.pipe( pE, fileOTUsPreTSV, c_fileProgFilterSamples, str(fileOTUsPreTSV).replace( c_strPre, "" ),
		[strFeatures] )[0]
	
	strBase = str(fileOTUsTSV).replace( sfle.c_strSufTSV, "" )
	strFileR = strBase + sfle.c_strSufR
	sfle.sop( pE, "ln -s", [[fileProgFileR], [True, strFileR]] )
	
	strFileTXT = strBase + sfle.c_strSufTXT
	afileCur = sfle.op( pE, c_fileProgMultivariateMetadata, [[True, strFileTXT], [fileOTUsTSV], "-m", strMetadatum, "-p"] )
	Depends( afileCur, strFileR )
	return ( afileCur + [fileOTUsTSV, strFileR] )

afileMaAsLin = maaslin( c_dirMaAsLin, c_fileOTUsPCL, fileDataPCL, c_strREFeatures, c_fileProgFileR, c_strMetadatum )
sfle.cmd( pE, c_fileProgMAssocs2Table, c_fileMAssocsTXT, [[afileMaAsLin[-1]], [afileMaAsLin[0]]] )
filePrePCL = sfle.pipe( pE, afileMaAsLin[-2], c_fileProgTranspose,
	sfle.d( fileDirTmp, sfle.rebase( afileMaAsLin[-2], True, c_strPre + sfle.c_strSufPCL ) ) )[0]
filePCL = Command( str(filePrePCL).replace( c_strPre, "" ), filePrePCL, sed( r'\t[AC]P', "\tP" ) )[0]
sfle.sink( pE, filePCL, c_fileProgPlotLevels,
	[[True, sfle.d( pE, fileDirOutput, sfle.rebase( c_fileMAssocsTXT, True, sfle.c_strSufPDF ) )],
	"-i", [c_fileMAssocsTXT], "-m", c_strMetadatum, "-k", c_strKeys] )

# #===============================================================================
# # RISK
# #===============================================================================

# c_strRERiskFeatures			= r'^(Exp.+)$'
# c_fileRiskMetadataPreTXT	= sfle.d( pE, fileDirTmp, sfle.rebase( c_fileInputRiskMetadataTXT, True, c_strPre + sfle.c_strSufTXT ) ) 
# c_fileRiskMetadataTXT		= str(c_fileRiskMetadataPreTXT).replace( c_strPre, "" ) 
# c_fileRiskOTUsPCL			= sfle.d( pE, fileDirTmp, sfle.rebase( c_fileInputRiskOTUsTXT, True, sfle.c_strSufPCL ) ) 
# c_dirRiskMaAsLin			= Dir( sfle.d( fileDirTmp, "risk" ) )

# sfle.pipe( pE, c_fileInputRiskMetadataTXT, c_fileProgTranspose, c_fileRiskMetadataPreTXT )
# sfle.cmd( pE, c_fileProgMergeTables, c_fileRiskMetadataTXT, [[c_fileRiskMetadataPreTXT], [c_fileInputRiskMapTXT]] )
# sfle.pipe( pE, c_fileInputRiskOTUsTXT, c_fileProgQiime2PCL, c_fileRiskOTUsPCL )
# apProcs = [
# 	sfle.CProcessor( "mtd",	None,	sfle.CCommand( c_fileProgMergeMetadata, [[c_fileRiskMetadataTXT], "-m 0"] ) ),
# 	sfle.CProcessor( "flt",	None,	sfle.CCommand( c_fileProgFilterFeatures, [c_strRiskMetadatum, c_dMean, c_dMax] ),
# 		sfle.CTarget( None, fileDirOutput ) ),
# ]
# fileRiskPCL = sfle.CProcessor.pipeline( pE, apProcs, c_fileRiskOTUsPCL )[0]

# #------------------------------------------------------------------------------ 

##maaslin( c_dirRiskMaAsLin, c_fileRiskOTUsPCL, fileRiskPCL, c_strRERiskFeatures, c_fileProgRiskFileR, c_strRiskMetadatum )
