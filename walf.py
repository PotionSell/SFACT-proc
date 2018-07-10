#!/usr/bin/python

import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from os.path import isfile, join, splitext, split, exists, sep
from os import listdir, makedirs


#packages for Pyraf support:
import sys
from io import StringIO
#from pyraf.iraf import splot, imhead
from pyraf import iraf
from time import sleep

def run1Dspec(specFile, objName):
    '''
    Processes a 1-dimensional .fits spectrum using ALFA.
    
    Args:
    specFile    absolute path to a 1 dimensional .fits spectrum file
    ????????????????????????????
    
    Returns:
    dataDict    dictionary with data keys: redshift; observedWavelength;
                physicalWavelength; lineFlux; equivalentWidth. Data entries are 
                either computed or derived using the ALFA program. Or, the dict 
                is empty if the user skips the spectrum.
                
    Output Files:
    /fittedspectra              directory holding .png plots of fitted spectra
    1dspec.*.fits.png           a .png Python plot of a fitted spectrum & its lines
    splot.log                   a log of line measurements made in IRAF's splot
                                (look up 'IRAF splot' for more details)
    /ALFAoutput                 directory holding ALFA's output files
    1dspec.*.fits_fit           (from ALFA) fitted spectrum and continuum data
    1dspec.*.fits_lines         (from ALFA) list of fitted lines
    1dspec.*.fits_lines.tex     (from ALFA) list of fitted lines (in LaTeX)
    '''

    specPath, specName = split(specFile)
    #get the true datapath, which is one directory up
    datapath = join(sep,*(specPath.split(sep)[:-1]))

    #directories for output files (create if nonexistent)
    outPath = join(datapath, 'ALFAoutput')
    imagepath = join(datapath, 'fittedspectra')
    if not exists(outPath):
        makedirs(outPath)
    if not exists(imagepath):
        makedirs(imagepath)
    
    
    #keep going until the user gets their desired good line measurement
    good = True
    dataDict = {}
    while(good):
        ##Obtain measurement of a single line to get spectrum's redshift.
        
        #run splot to let user measure 1 line
        logFile = join(specPath, 'splot.log')
        print('\n*****In the splot window, measure one line with a known rest ' \
                'wavelength, then press q. If there are no lines, just press q.')
        iraf.splot( images=specFile, save_file=logFile )

        #prompt user for the known reference wavelength. Otherwise, skip.
        while(True):
            try:
                #provide shortcuts for certain lines, using high precision
                key = eval(input(
                                '*****If you measured a line, enter:\n'             \
                                '0 for no lines - skip this spectrum\n'             \
                                '1 for H-alpha \n'                                  \
                                '2 for OIII \n'                                     \
                                '3 for OII \n'                                      \
                                'the wavelength of another line\n'                  \
                                ))
                key = float(key)        #trigger a Name/ValueError if not a float
            except (NameError, ValueError, SyntaxError):
                print('Invalid input - please enter 0, 1, 2, 3, or a wavelength:\n')
                continue
                
            if key == 0: return {}
            elif key == 1: restW = 6562.82      #H-alpha
            elif key == 2: restW = 5006.84      #OIII
            elif key == 3: restW = 3727.53      #OII blend
            else: restW = key
            break


        #read splot.log to get the wavelength measured by user in splot
        cols = ['center','cont','flux','eqw','core','gfwhm','lfwhm']
        splotData = pd.read_table( logFile, names=cols, keep_default_na=False, sep='\s+' )
        #mark the invalid rows...
        splotData = splotData.apply(pd.to_numeric, errors='coerce')
        #...and get rid of them
        splotData = splotData.dropna(subset = ['center'])
        #get the most recent addition to the log file
        obsW = splotData.iloc[-1][0]

        z = obsW / restW - 1.
    
    
        ##Call ALFA to fit the spectrum and lines.
        alfapath = '/home/bscousin/software/bin/alfa'
        callALFA(alfapath, z, specFile, outPath)
        
        ##Extract info from ALFA's output for plotting and calculations.

        #read in spectral data from ALFA
        fitFname = join(outPath, specName+'_fit')
        specData = np.genfromtxt(fitFname, skip_header=2)
        specW = specData[:,0]
        specF = specData[:,1]
        modelF = specData[:,2]
        contF = specData[:,4]

        #read in fitted-line data from ALFA
        linesFname = join(outPath, specName+'_lines')
        mycols = ['obsW', 'realW', 'flux', 'uncertainty', 'peak', 'fwhm']
        rawLines = pd.read_csv(linesFname, sep='\s+', names=mycols)
        #convert to dataframe of floats, ignoring strings
        rawLines = rawLines.apply(pd.to_numeric, errors='coerce')
        #find lines that are actual detections, not just upper limits
        rowidx = np.where(~np.isnan(rawLines.uncertainty))
        lineData = rawLines.iloc[rowidx]

        nLines = np.shape(lineData)[0]
        obsW = lineData['obsW']
        realW = lineData['realW']
        flux = lineData['flux']
        uncertainty = lineData['uncertainty']
        peak = lineData['peak']
        fwhm = lineData['fwhm']
        
        
        ##Compute equivalent widths.
        
        #get indices of lines, and find the average continuum around each
        disp = 1.4076
        numFWHM = 2.5
        idx = np.empty(nLines, dtype=int)
        contAvg = np.empty(nLines, dtype=np.float64)    #flux values are very small, so need high precision!

        def firstIndex(a, val, tol):
            '''Gets the first index in a numpy array (a) that's within a 
            tolerance (tol) of some value (val)'''
            idx = (np.abs(a - val) <= tol).argmax()
            return idx

        for i in range(nLines):
            idx[i] = firstIndex(specW, obsW.iloc[i], disp)
            #idx is close (+-1) to the exact index - this is okay since we only use it for 
            #getting averages across a range of indices anyways
            
            idxRange = int(np.rint( fwhm.iloc[i] * numFWHM / disp / 2))
            idxWindow = [idx[i]-idxRange, idx[i]+idxRange]
            
            contAvg[i] = np.mean( contF[ idxWindow[0] : idxWindow[1] ] )

        eqWidth = np.array(peak / contAvg)


        ##Create spectrum plot.
        
        plt.figure(figsize=(18,12))
        plt.plot(specW, specF, c='grey', label='original spectrum')
        plt.plot(specW, contF, c='black', ls='--', label='continuum fit')

        ax = plt.gca()
        ymin, ymax = ax.get_ybound()
        for i in range(nLines):
            plt.axvline(x=obsW.iloc[i], c='lightgrey', ls='--')

        #label peak line fluxes, individually with wavelength info
        plt.plot(obsW, peak + contAvg, 'o', c='orchid', label='estimated line peak')
        for i, j, k in zip(obsW, peak+contAvg, realW):
            #shift the labels a bit by a correction so they aren't on top of 
            #data points. Shift is intelligently scaled to the axis size.
            xcorr = max(obsW)/-100.
            ycorr = max(peak+contAvg)/30.
            ax.annotate(str(i), xy=(i+xcorr,j+ycorr))
            ax.annotate('('+str(k)+')', xy=(i+xcorr,j+ycorr/2.75))

        plt.title('Fitted lines for spectrum ' +specName+ '; object ID ' +objName, size=24)
        plt.xlim(min(specW-100), max(specW+100))
        plt.xlabel(r'wavelength ($\AA$)', size=20)
        plt.ylabel(r'flux ($units$)', size=20)
        plt.legend(loc='upper left', prop={'size': 16})
        
        plt.savefig( join(imagepath, specName+'.png') )
        print('The fitted spectrum has been saved in \'fittedspectra.\'\n')
        plt.ion()
        plt.show()
#        import pdb; pdb.set_trace()

        print('Redshift of measured line: ', round(z,4))
        
        #Check if user wants to redo the line selection and ALFA call.
        #if they redo, then break to the main loop and try again
        #if they skip, then break to the main loop and do not try again, returning empty
        #if they keep, then break to the main loop and do not try again, returning data
        while(True):     #loop until user gives valid input
            try: 
                redo = eval(input(
                            '*****Would you like to redo the line estimation? ' \
                            'Enter 0 to skip; 1 to redo; 2 to keep the current estimate:\n'
                            ))
                redo = int(redo)        #trigger a ValueError if user gives bad input
            except (ValueError, NameError, SyntaxError):     #user gives non-integer input
                print('Invalid input - please enter 0, 1, or 2:\n')
                continue

            if (redo != 0 and redo != 1 and redo != 2): 
                print('Invalid input - please enter 0, 1, or 2:\n')
                continue
            elif redo == 0:
                good = False
                #set the return value to be empty
                dataDict = {}
                return dataDict
            elif redo == 2:
                good = False      
                ##set the return value to our data
                dataDict = {
                            'redshift': z, 'observedWavelength': np.array(obsW),\
                            'physicalWavelength': np.array(realW),              \
                            'lineFlux': np.array(flux),                         \
                            'uncertainty': np.array(uncertainty),               \
                            'peak': np.array(peak),                             \
                            'fwhm': np.array(fwhm),                             \
                            'equivalentWidth': eqWidth
                            }
                return dataDict
            #and then break to the main loop
            break
    return dataDict


def runMultispec(fpath):
    '''
    Identifies and measures the wavelengths, fluxes, and equivalent widths of 
    emission lines in a 2-dimensional set of spectra using Pyraf and ALFA.
    
    Prompts the user to measure and identify one line, via Pyraf's/Iraf's splot,
    in order to derive a redshift to be used in ALFA.
    
    Args:
    fpath       absolute path to a 2-dimensional set of spectra .fits file
    
    Returns:
    -----
    
    Output Files:
    /1dspectra      directory for 1D spectra generated from the input 2D spectrum
    1dspec.*.fits   a single 1D spectrum generated from the input 2D spectrum; 
    *.txt           output file containing fields for all spectra that the user
                    does not skip
    '''

    datapath, basename = split(fpath)
    fname, __ = splitext(basename)
    
    #directories for output files (create if nonexistent)
    specPath = join(datapath, '1dspectra')
    if not exists(specPath):
        makedirs(specPath)
    
    ##Prepare to write to file. Each spectrum is processed and written with the 
    ##output file open, so processed data is still written even if the code 
    ##crashes/escapes. 
    with open( join(datapath, fname)+'.txt', 'w+' ) as f:
        ##Access 2Dspectrum and begin processing.

        #hacky workaround to get output printed by IRAF's imhead as a string
        old_stdout = sys.stdout
        sys.stdout = mystdout = StringIO()
        iraf.imhead(fpath, longheader='yes')
        sys.stdout = old_stdout

        #access strings from stdout, and format as a numpy array
        imheadStr = mystdout.getvalue()
        imheadStr = imheadStr.split('\n')
        imheadStr = np.array(imheadStr)

        #find rows for APID and APNUM (yes, they're different rows)
        apidx1 = np.where(['APID' in x for x in imheadStr])
        apidx2 = np.where(['APNUM' in x for x in imheadStr])

        apIDRows = imheadStr[apidx1]
        apNumRows = imheadStr[apidx2]

        #split of the different entries in the row...
        apNumRowValues = [x.split('\'') for x in apNumRows]
        #...and select the part of the row that contains the # and flag...
        nums = [x[1] for x in apNumRowValues]
        #...and select only the part containing the #...
        apNums = np.array( [ x.split()[0] for x in nums ] )
        #and the part containing the flag
        apFlags = np.array( [ int( x.split()[1] ) for x in nums ], dtype=bool )
        
        #also extract object names/IDs
        apIDRowValues = [x.split('\'') for x in apIDRows]
        names = [x[1] for x in apIDRowValues]
        apNames = np.array( [x.split()[0] for x in names] )
        
        goodNums = apNums[apFlags]
        goodNames = apNames[apFlags]
        
        nSpec = len(goodNums)
#        import pdb; pdb.set_trace()
        
        ##Start writing/processing, spectrum by spectrum:
        f.write('ID\t\tredshift\tobservedW\tphysicalW\tflux\t\tuncertainty\tpeak\t\tfwhm\teqWidth\n')   
        for i in range(nSpec):
#        for i in range(33, nSpec):
            ##Run scopy to make a file for the current spectrum.
            objName = goodNames[i]
            curAp = int(goodNums[i])
            pad_ap = format(curAp, '04d')
            specName = '1dspec.'+pad_ap+'.fits'
            specFile = join(specPath, specName)
            
            iraf.scopy(input=fpath, output=join(specPath, '1dspec'), apertures=curAp)
            
            ##Process the 1D spectrum
            data = run1Dspec(specFile, objName)
            
            #if user skipped the spectrum, data is empty, so continue to next spectrum
            if not data: continue
                        
            ##Proceed with writing the output file.
            print('Writing line data for: ' +specFile+ '\n')
            nLines = len(data['observedWavelength'])
            for j in range(nLines):
                f.write(objName                                                 \
                        +'\t' +str(round(data['redshift'], 5))                  \
                        +'\t\t'+ str(data['observedWavelength'][j])             \
                        +'\t\t'+ str(data['physicalWavelength'][j])             \
                        +'\t\t'+ str('%.2E' %data['lineFlux'][j])               \
                        +'\t'+ str('%.2E' %data['uncertainty'][j])              \
                        +'\t'+ str('%.2E' %data['peak'][j])                     \
                        +'\t'+ str(data['fwhm'][j])                             \
                        +'\t'+ str(np.round( data['equivalentWidth'][j], 2 ))   \
                        +'\n' 
                        )
        f.write('\n')
    
#    #scrape up all the files written by ALFA
#    allFiles = [f for f in listdir( join(datapath,savefolder) ) if f.endswith('.fits_fit') or f.endswith('.fits_lines')]


def run1Dfrom2D(fpath, apNum):
    '''
    Run scopy on a spectrum specified by apNum, and process it with run1Dspec()
    '''
    pass


def callALFA(alfapath, z, specFile, outPath):
    '''
    Internal function for run1Dspec. Calls the Automated-Line-Fitting-Algorithm
    (ALFA) from the command line using determined parameters.
    
    Args:
    alfapath    absolute path to the system's ALFA install (***may not be needed
                for proper installs***)
    z           redshift of the spectrum
    specFile    absolute path to the spectrum file
    outPath     absolute path to the output directory
    
    Returns:
    -----
    '''
    
    v = z * 299792	#km/s
    strongCat = ' optical_strongCustom.cat '
    deepCat = ' optical_deepCustom2.cat '
    #set window size used during continuum fitting
    try:
       winSize=eval(input('Enter the window size to be used for continuum fitting '\
                    '(an odd # of pixels - leave empty for default 101 pixels)"\n'))
    except SyntaxError:
        winSize = 100

    paramStr =  ' -ul  -n 0  -vg ' +str(v)+ ' -cw ' +str(winSize)+ \
                ' -o ' +outPath+ ' -sc' +strongCat+ '-dc' +deepCat+ ' '
    #-ul gets upper limit output; -n specifies normalization (none in this case)
    # -vg specifies velocity guess; -cw specifies the window size
    #-o specifies the path for output files; -sc/-dc specify line catalogs
    #See ALFA's manual for more info: >>> man alfa

    alfaCommand = alfapath + paramStr + specFile
#    process = subprocess.call(alfaCommand, shell=True)
    process = subprocess.call(alfaCommand, shell=True, stdout=subprocess.PIPE)	#muted output




fpath = eval("input('Enter the full path of 2D spectrum fits file: ')")
fpath = '/home/bscousin/iraf/Team_SFACT/hadot055A/hadot055A_comb_fin.ms.fits'
#fpath = '/home/bscousin/iraf/Team_SFACT/hadot055A/hadot055A_comb_fp.ms.fits'
runMultispec(fpath)


