#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 23 11:43:07 2022

@author: dan
"""

import os
import sys

sys.path.append('wma_pyTools/')
startDir=os.getcwd()
#some how set a path to wma pyTools repo directory
#wmaToolsDir='wma_pyTools'
#wmaToolsDir='..'
#os.chdir(wmaToolsDir)
print(os.getcwd())
print(os.listdir())
import wmaPyTools.roiTools
import wmaPyTools.analysisTools
import wmaPyTools.segmentationTools
import wmaPyTools.streamlineTools
import wmaPyTools.visTools

#os.chdir(startDir)

import os
import json
import numpy as np
import nibabel as nib
import pandas as pd

import wmaPyTools.genUtils
if wmaPyTools.genUtils.is_docker():
    import matplotlib 
    #trying workaround described here:
    #https://github.com/matplotlib/matplotlib/issues/18022
    #import matplotlib.backends
    matplotlib.use('Agg')


# load inputs from config.json
with open('config.json') as config_json:
	config = json.load(config_json)
    
inflateParam=config['inflateParam']

outDir='output'
if not os.path.exists(outDir):
    os.makedirs(outDir)
if not os.path.exists(os.path.join(outDir,'images')):
    os.makedirs(os.path.join(outDir,'images'))
#set to freesurfer output path for this subject
fsPath=config['freesurfer']
#get the parc input, if its there
parcIn=config['parc']

#you may need to convert the .mgz files to .nii.gz using the mr_convert command
#also, you may need to rename the subsequent aparcDk atlas file to it's standard name:

#read in the input parc, from the relevant source, if available
if not np.logical_or(parcIn=='',parcIn==None):    
    inputAtlas=nib.load(parcIn)
    #do a bit of preprocessing
    inputAtlas=wmaPyTools.roiTools.preProcParc(inputAtlas,deIslandBool=True,inflateIter=inflateParam,retainOrigBorders=False,maintainIslandsLabels=None,erodeLabels=None)    
    
    lookupTable=wmaPyTools.genUtils.parcJSON_to_LUT(config['label'])
    
elif  not np.logical_or(fsPath=='',fsPath==None):
    atlasName='aparc.a2009s+aseg'
    try:
        inputAtlas=nib.load(os.path.join(fsPath,'mri/'+atlasName+'.nii.gz'))
    except:
        #can nibael handle mgz?
        inputAtlas=nib.load(os.path.join(fsPath,'mri/'+atlasName+'.mgz'))
    #in either case get the lookup table
    lookupTable=pd.read_csv('FreesurferLookup.csv')
    #and do preprocessing
    inputAtlas=wmaPyTools.roiTools.preProcParc(inputAtlas,deIslandBool=True,inflateIter=inflateParam,retainOrigBorders=False,maintainIslandsLabels=None,erodeLabels=[2,41])    
    
else:
    inputAtlas=None
    
refAnatT1=config['anat']
if not np.logical_or(refAnatT1=='',refAnatT1==None):
    refAnatT1=nib.load(refAnatT1)
    #robust load; sometimes hi res T1s have very strange affines
    refAnatT1 = nib.nifti1.Nifti1Image(refAnatT1.get_data(), np.round(refAnatT1.affine,4), refAnatT1.header)
else:
    #set refAnatT1 to none
    refAnatT1=None

#set to path to target whole brain tractogram
#smaller = faster
tractogramPath=config['tractogram']
tractogramLoad=nib.streamlines.load(tractogramPath)
#is this creating inf values?
#streamlines=wmaPyTools.streamlineTools.orientAllStreamlines(tractogramLoad.streamlines)
streamlines=tractogramLoad.streamlines

#load the wmc
classification=wmaPyTools.streamlineTools.matWMC2dict(config['wmc'])
#reminder
#wmc_Dict['names']=tractNames
#wmc_Dict['index']=indices.tolist()

outJsonDict={'images':[]}
for tractIterator,iTractName in enumerate(classification['names']):
    
    print ('plotting figures for ' + iTractName)
    #get the current name
    currentName=iTractName
    #get the bool vec of the streamline indexes
    currentIndexesBool=[iIndex==tractIterator+1 for iIndex in classification['index']]
    #check to make sure something is there
    if not np.sum(currentIndexesBool)==0:
        #make an output directory for this tract
        currFigOutDir=os.path.join(outDir,'images',currentName)
        if not os.path.exists(currFigOutDir):
            os.makedirs(currFigOutDir)
        #create the requested visualizations
        wmaPyTools.visTools.multiPlotsForTract(streamlines[currentIndexesBool],atlas=inputAtlas,atlasLookupTable=lookupTable,refAnatT1=refAnatT1,outdir=currFigOutDir,tractName=currentName,makeGifs=config['gifFlag'],makeTiles=config['tileFlag'],makeFingerprints=config['fingerprintFlag'],makeSpagetti=config['spagettiFlag'])
        #generate a json info dict for the requested images
        currentTractDict=wmaPyTools.visTools.jsonFor_multiPlotsForTract(saveDir=currFigOutDir,tractName=currentName,makeGifs=config['gifFlag'],makeTiles=config['tileFlag'],makeFingerprints=config['fingerprintFlag'],makeSpagetti=config['spagettiFlag'])
        #append it to the dictionary
        outJsonDict['images']=outJsonDict['images']+currentTractDict['images']
    
with open(os.path.join(outDir,"images.json"), "w") as outfile:
    #dump or dumps?  I have no idea
    json.dump(outJsonDict, outfile)

print ('figure generation complete')
