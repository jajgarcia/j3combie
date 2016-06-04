#!/usr/bin/python
#
# j3combine.py
#
# Combination procedure designed by J. Garcia, J. Steiner,
# and J. McClintock (see Garcia et al. 2015, ApJ 813, 84G)
# 
# Observatinons are chosen from the given input file
# (must be listed one per line). All data is loaded and
# the model is fit using PyXspec. The standard procedure uses
# a simple absorbed power-law model (Tbabs*pow), but other 
# models are possible. The individual residuals are combined,
# and a fake power-law continuum is added based on the average
# parameters from all the infividual fits.
#
#  TO-DO:
# - Include an option to define the obspath
# - Make the code to create the logs folder if it doesn't exist
# - Also, include an option to turn on/off the logs
#
# Requires: xspec, astropy, time
#
import sys 
from xspec import *
from optparse import OptionParser
import os,os.path
import glob
from astropy.io import fits as pyfits
from astropy.time import Time
import time
from datetime import datetime
# ------------------------------------------------------------------------------
# regrid data set
#
# mesh - original mesh
# dat - original data
# pmesh - new mesh
#
# return - array of data on new mesh
def regrid_dat(mesh,dat,pmesh):
  out=[]
  for ix in range(len(pmesh)):
    x=pmesh[ix]

    if x < mesh[0]:            # interpolate to the left
      out.append(interpolate(x,mesh[0],dat[0],mesh[1],dat[1]))
      #out.append(0.1)
    elif x > mesh[-1]:         # interpolate to the right
      out.append(interpolate(x,mesh[-2],dat[-2],mesh[-1],dat[-1]))
      #out.append(0.1)
    else:                      # normal interpolation
      ixp=0
      while mesh[ixp] < x: ixp+=1
      out.append(interpolate(x,mesh[ixp-1],dat[ixp-1],mesh[ixp],dat[ixp]))

  return out
# ------------------------------------------------------------------------------
#
# perform linear interpolation (or extrapolation)
#
# x0 - where to interpolate
# x1 - 1st x point
# y1 - 1st y point
# x2 - 2nd x point
# y2 - 2nd y point
#
# return - interpolated y point
def interpolate(x0,x1,y1,x2,y2):
  if (x1 == x2):
    print "cannot interpolate given indentical x1 and x2"
    raise
  slope=(y2-y1)/(x2-x1)
  return y1+slope*(x0-x1)
#
# ------------------------------------------------------------------------------
#
# Read a file and return a list of lines
#
# infile - File name
# lines  - List of lines in the file
#
def read_lines(infile):
  lines=[]
  fin=open(infile,'r')
  while 1:
    tline=fin.readline()
    if tline == '': break                # Reaches EOF
    if tline.strip() == '': continue     # Avoid empty lines
    lines.append(tline.strip().split())
  fin.close
  return lines
#
# ------------------------------------------------------------------------------
#
# MAIN PROGRAM
#
#
#
version='0.1a'
date='- Wed May 25 18:59:34 EDT 2016 -'
author='Javier Garcia'
#
ul=[]
ul.append("usage: %prog [options] PREFIX")
ul.append("")
ul.append("Get total counts in different bands for a given observation")
ul.append("PREFIX can be a single PHA file or a group (e.g. *.pha)")
usage=""
for u in ul: usage+=u+'\n'

parser=OptionParser(usage=usage)
parser.add_option("-v","--version",action="store_true",dest="version",default=False,help="show version number")
parser.add_option("-p","--prefix",dest="prefix",default="spec",help="specify prefix for output file")
parser.add_option("-o","--output",dest="outfile",default="",help="specify alternative output file")

(options,args)=parser.parse_args()

if options.version:
  print 'j3combine.py version:',version,date
  print 'Author:',author
  sys.exit()

if len(args) == 0:
  parser.print_help()
  sys.exit(0)

prefix=options.prefix
finalsp=options.outfile
if finalsp == "":
  finalsp=prefix+'_J3combined.pha'

if os.path.isfile(finalsp):
  print 'Error: output file',finalsp,'already exist!'
  print 'Aborting...'
  sys.exit()

# Get current universal date and time 
currtime=str(datetime.utcnow())

#-----

# No chatter
Xset.chatter = 0
#Xset.logChatter = 0

# Abundances and cross sections
Xset.abund = "wilm"
Xset.xsect = "vern"

# Query
Fit.query = 'yes'

# Load local models
#AllModels.lmod("jscrab")

# Get current directory
currpath = os.getcwd()

# Observations path
obspath='/Users/javier/H1743-322/observations/'

# List of spectrum files
files = read_lines(args[0])

# Change dir to observations path
os.chdir(obspath)

#----- LOOP OVER OBSERVATIONS ---#
residuals=[]
data=[]
background=[]
refmesh=False
gammas=[]
norms=[]
sotimes=[]
batimes=[]
for tmpfile in files:

  specfile = tmpfile[0]
  print specfile

  # Change dir to observations path
  os.chdir(obspath)

  # Check if specfile exist
  if not os.path.isfile(specfile):
    print 'Warning: spectrum file',specfile,'does not exist!'
    print 'Skiping...'

  else:  # 

    # 
    # Read spectrum file
    hdulist = pyfits.open(specfile)

    # Get times
    dateobs=hdulist[1].header['DATE-OBS']
    time = Time(dateobs, format='isot', scale='utc')
    timemjd=time.mjd

    # Get response file name
    respfile=hdulist[1].header['respfile']
    if not os.path.isfile(respfile):
      print 'Error: response file',specfile,'does not exist!'
      print 'Aborting...'
      sys.exit()
    # Here I need to set the response file for the combined spectrum

    #Background
    backfile=hdulist[1].header['backfile']
    if not os.path.isfile(backfile):
      print 'Error: response file',backfile,'does not exist!'
      print 'Aborting...'
      sys.exit()

    # Get the exposure for background
    backlist = pyfits.open(backfile)
    soet = hdulist[1].header['exposure']   # Source Exposure Time
    baet = backlist[1].header['exposure']  # Background Exposure Time
 
    # Save the times for combination later
    sotimes.append(soet)   
    batimes.append(baet)

    # Load data
    AllData('1:1 '+specfile)
    s1 = AllData(1)

    # Go back to the working directory
    os.chdir(currpath)

    # Ignore/notice data
    s1.ignore("1-4,45.-**")

    # Define the Model
    AllModels += "const*tbabs*pow"
    m1 = AllModels(1)

    m1(1).values = "1 -1"    # Constant
    m1(2).values = "0.7 -1"  # nH
    m1(3).values = "1.5 1"   # Gamma
    m1(4).values = "1. 1"    # Normalization

    # Fit
    Fit.renorm()
    Fit.perform()

    # Create and open a log file
    logFile = Xset.openLog('logs/fit-'+specfile+'.log')
    logFile = Xset.log

    # Calculate Errors
    Fit.error("maximum 200. 2.706 3 4")

    # Close currently opened log file.
    Xset.closeLog()

    # Equivalent to show all
    logFile = Xset.openLog('logs/fit-'+specfile+'.txt')
    logFile = Xset.log
    s1.show()
    m1.show()
    Fit.show()
    Xset.closeLog()
  
    # Save Gamma and Normalization (and their uncertainties)
    gammas.append([m1(3).values[0],m1(3).error[0],m1(3).error[1]])
    norms.append([m1(4).values[0],m1(4).error[0],m1(4).error[1]])

    #
    s1.notice("**")
    Plot.xAxis = "keV"
    Plot.background = True
    Plot("Counts")
    dvals = Plot.y()            # Data: Data-Model (counts/keV)
    bkg = Plot.backgroundVals() # Backgroud: (counts/keV)
    Plot("residuals")
    xvals = Plot.x()            # Center bin energy (keV)
    yvals = Plot.y()            # Residuals: Data-Model (counts/sec/keV)
    xErrs = Plot.xErr()         # Half of bin width (keV)
    yErrs = Plot.yErr()         # Sigma: Error bar (counts/sec/keV)

    # Save the reference grid, data, residuals, and background
    if refmesh is False:
       energies=[]
    ores=[]
    odata=[]
    oback=[]
    for i in range(len(xvals)):
      if refmesh is False:
         energies.append(xvals[i])             # Reference energy grid (just choose the first one)
      oback.append(bkg[i]*2.*xErrs[i])         # Background (counts)
      ores.append(yvals[i]*soet*2.*xErrs[i])   # Residuals (counts)
      odata.append(dvals[i]*2.*xErrs[i])       # Data (counts)

    # Map the observation to the reference grid
    nres = regrid_dat(xvals,ores,energies)
    ndata = regrid_dat(xvals,odata,energies)
    nback = regrid_dat(xvals,oback,energies)

    # I am not sure if I need to save the individual spectra
    # (might take lots of RAM!) but just in case for now
    residuals.append(nres)
    data.append(ndata)
    background.append(nback)

    # Define the reference response and background
    # (just take the first one in the list)
    if refmesh is False:
       refresponse=respfile
       refsource=specfile
       refmesh=True
    # Unload data
    AllData -= s1

    # Unload the model
    AllModels.clear()
#-----

# Combine residuals from all observations
cres = [sum(x) for x in zip(*residuals)]
cdata = [sum(x) for x in zip(*data)]
cback = [sum(x) for x in zip(*background)]

# Output: ratio spectrum
csig=[]
for i in range(len(energies)):
   csig.append((cdata[i]+(soet/baet)*cback[i]+(soet/baet)**2.*cback[i])**0.5)

# Average Gamma, Normalization and times (to be used to create fake continuum)
avgGamma = [sum(x)/len(gammas) for x in zip(*gammas)]
avgNorm = [sum(x)/len(norms) for x in zip(*norms)]
totsource = sum(sotimes)
totback = sum(batimes)

#### Simulate the fake continuum ####

# Change dir to observations path
os.chdir(obspath)

# Define the Model
m1 = Model("tbabs*pow")

m1(1).values = "0.7 -1"               # nH
m1(2).values = str(avgGamma[0])+" 1"  # Gamma
m1(3).values = str(avgNorm[0])+" 1"   # Normalization

#response, arf, background, exposure, correction, backExposure, fileName
#fs1 = FakeitSettings(str(refresponse),"",str(refbackground),str(totsource),1.,str(totback),"fake-pca.pha")
# Not sure if I need to simulate the backgorund!!!
fs1 = FakeitSettings(refresponse,"","",totsource,1.,0.,"fake-pca.pha")
AllData.fakeit(1, fs1)

s1 = AllData(1)
Plot.xAxis = "keV"
Plot("Counts")
xErrs = Plot.xErr()  # Half of bin width (keV)
plvals = Plot.y()    # Power-law (counts/keV)

# Unload data
AllData.clear()

# Go back to the working directory
os.chdir(currpath)

# Add the combined residuals and the fake power-law continuum
spcomb = [x + y*2.*z for x, y, z in zip(cres, plvals, xErrs)]

# Make sure the final spectrum counts are all integers
spcomb = map(int, spcomb)

# Output
#
#
tmpspec = pyfits.open(obspath+refsource)
tmpspec_hdu = tmpspec['SPECTRUM']
#
tmpspec_hdu.header["EXPOSURE"] = totsource
tmpspec_hdu.header["BACKFILE"] = 'NONE'
tmpspec_hdu.header["CORRFILE"] = 'NONE'
#
tmpspec_counts = tmpspec_hdu.data.field("COUNTS")
tmpspec_staterr = tmpspec_hdu.data.field("STAT_ERR")
tmpspec_counts[0:128] = spcomb[0:128]
tmpspec_staterr[0:128] = csig[0:128]
#
# Include HISTORY comments
header = tmpspec[0].header      # reads header
header['HISTORY'] = 'Observation file created by j3combine.py version '+version
header['HISTORY'] = 'Files combined: '+str(files)
header['HISTORY'] = 'Creation date: '+currtime+' UTC'
#
tmpspec.writeto(finalsp)
#
sys.exit()
# ------------------------------------------------------------------------------
