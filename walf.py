#!/usr/bin/python

import matplotlib
matplotlib.use('TkAgg')     #graphics setting so that iraf.splot doesn't blow up on Macs

import subprocess
import collections
import math
from os.path import isfile, join, splitext, split, exists, sep
from os import listdir, makedirs
from time import sleep

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#packages for Pyraf support:
import sys
from io import StringIO
from pyraf import iraf


def run1Dspec(specFile, objID, deltaW=0):
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
    winSize = 150
    good = True
    while(good):
        userW, restW, skip = promptLine(specFile, 'science')
        
        if skip:
            objFlag = classifySpectrum()
            objDF = sparseObjDF(objID, objFlag)
            lineDF = sparseLineDF()
            return objDF, lineDF
        
        #adjust for sky line correction
        userW = userW - deltaW
        
        z = userW / restW - 1.
        
        print('\nWavelength of measured line: ', userW)
        print('Redshift of measured line: ', round(z,5), '\n')

        ##Call ALFA to fit the spectrum and lines.
        callALFA(z, winSize, specFile, outPath)
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
        mycols = ['observeW', 'physW', 'flux', 'sigmaFlux', 'peak', 'fwhm']
        rawLines = pd.read_csv(linesFname, sep='\s+', names=mycols)
        #convert to dataframe of floats, ignoring strings
        rawLines = rawLines.apply(pd.to_numeric, errors='coerce')
        #find lines that are actual detections, not just upper limits
        rowidx = np.where(~np.isnan(rawLines.sigmaFlux))
        goodLines = rawLines.iloc[rowidx]

        nLines = np.shape(goodLines)[0]
        observeW = goodLines['observeW']
        physW = goodLines['physW']
        flux = goodLines['flux']
        sigmaFlux = goodLines['sigmaFlux']
        peak = goodLines['peak']
        fwhm = goodLines['fwhm']
        
        
        ##Compute equivalent widths.
        
        #get indices of lines in the spectrum vector, and find the average 
        #continuum around each
        disp = 1.4076
        numFWHM = 2.5
        idx = np.empty(nLines, dtype=int)
        contAvg = np.empty(nLines, dtype=np.float64)    #flux values are very small, so need high precision!
        
        matchIdx = lambda a,val,tol : np.abs(a-val) <= tol
#        def firstIndex(a, val, tol):
#           '''Gets the first index in a numpy array (a) that's within a 
#            tolerance (tol) of some value (val)'''
#            idx = (np.abs(a - val) <= tol).argmax()
#            return idx

        for i in range(nLines):
            idx[i] = matchIdx(specW, observeW.iloc[i], disp).argmax()
#            idx[i] = firstIndex(specW, observeW.iloc[i], disp)
            #idx is close (+-1) to the exact index - this is okay since we only use it for 
            #getting averages across a range of indices anyways
            
            idxRange = int(np.rint( fwhm.iloc[i] * numFWHM / disp / 2))
            idxWindow = [idx[i]-idxRange, idx[i]+idxRange]
            
            contAvg[i] = np.mean( contF[ idxWindow[0] : idxWindow[1] ] )
        eqWidth = np.array(flux / contAvg)


        ##Create spectrum plot.
        
        #create plot size depending on user's monitor
        #window = plt.get_current_fig_manager().window
        #screenX, screenY = window.wm_maxsize()
        #dpi = 100
        #fig = plt.figure( figsize=( int(screenX/dpi/1.5), int(screenY/dpi/1.2)),dpi=dpi )
        fig = plt.figure(figsize=(14,8))
        
        #manually add toolbar
        #from matplotlib.backends.backend_tkagg import NavigationToolbar2TkAgg as NavigationToolbar
        #win = fig.canvas.manager.window
        #canvas = fig.canvas
        #toolbar = NavigationToolbar(canvas, win)

        plt.plot(specW, specF, c='grey', label='original spectrum')
        plt.plot(specW, contF, c='black', ls='--', label='continuum fit')

        ax = plt.gca()
        ymin, ymax = ax.get_ybound()
        for i in range(nLines):
            plt.axvline(x=observeW.iloc[i], c='lightgrey', ls='--')

        #label peak line fluxes, individually with wavelength info
        plt.plot(observeW, peak + contAvg, 'o', c='orchid', label='estimated line peak')
        for i, j, k in zip(observeW, peak+contAvg, physW):
            #shift the labels a bit by a correction so they aren't on top of 
            #data points. Shift is intelligently scaled to the axis size.
            xcorr = max(observeW)/-100.
            ycorr = max(peak+contAvg)/30.
            text = str(i) +'\n('+ str(k) +')'
            ax.annotate(text, xy=(i+xcorr,j+ycorr))
#        import pdb; pdb.set_trace()
        plt.title('Spectrum ' +specName+ '; object ID ' +objID, size=24)
        plt.xlim(min(specW-100), max(specW+100))
        plt.xlabel(r'wavelength ($\AA$)', size=20)
        plt.ylabel('flux', size=20)
        plt.legend(loc='upper left', prop={'size': 16})
        
        plt.savefig( join(imagepath, specName+'.png') )
        print('The fitted spectrum has been saved in \'fittedspectra.\'\n')
        
        import pickle; pickle.dump(fig, open( join(imagepath,specName+'.pickle'), 'wb'))
        plt.ion()
        plt.show()
        
        #Check if user wants to redo the line selection and ALFA call.
        #if they redo, then break to the main loop and try again
        #if they skip, then break to the main loop and do not try again, returning basic global data
        #if they keep, then break to the main loop and do not try again
        while(True):     #loop until user gives valid input
            try: 
                validRedo = [0,1,2,3]
                redo = eval(input(
                            '*****The spectrum has been analyzed. Enter:\n'
                            '0 to skip (flag is saved, but not redshift or line data)\n'
                            '1 to skip (flag and redshift are saved, but not line data)\n'
                            '2 to remeasure the spectrum\n'
                            '3 to keep the current estimate (line data is saved)\n'
                            ))
                redo = int(redo)        #trigger a ValueError if user gives bad input
            except (ValueError, NameError, SyntaxError):     #user gives non-integer input
                print('Invalid input - please enter 0, 1, 2, or 3:\n')
                continue
            if redo not in validRedo: 
                print('Invalid input - please enter 0, 1, 2, or 3:\n')
                continue
            
            elif redo == 0:
                #return, ending the function
                objFlag = classifySpectrum()
                objDF = sparseObjDF(objID, objFlag)
                lineDF = sparseLineDF()
                return objDF, lineDF
            elif redo == 1:
                #return, ending the function
                objFlag = classifySpectrum()
                objDF = sparseObjDF(objID, objFlag, z)
                lineDF = sparseLineDF()
                return objDF, lineDF
            elif redo == 2:
                good = True
                #set window size used during continuum fitting
                while(True):
                    try:
                        winSize=eval(input('Enter the window size to be used for continuum fitting, '\
                                    'or press enter for the default of 150 pixels"\n'))
                        winSize = int(winSize)    #trigger a Name/ValueError if not an int
                    except (NameError, ValueError):
                        print('Invalid input - please enter an integer window size:\n')
                        continue
                    except SyntaxError:
                        winSize = 150
                    if winSize <= 0:
                        print('Invalid input - please enter an integer window size:\n')
                        continue
                    break
                
                break
            elif redo == 3:
                good = False
                objFlag = classifySpectrum()
                break

    ##Prepare output.
    nLines = len(observeW)
    #apply sky line correction to observed wavelengths
    observeW = observeW - deltaW
    #derive individual line redshifts
    lineZ = observeW/physW - 1.
    
    #compute an average redshift from data above a threshold SN ratio
    SNr = flux/sigmaFlux
    SNthres = 5.
    goodIdx = np.where(SNr >= SNthres)
    bestZ = np.mean(lineZ.iloc[goodIdx])
    bestZstd = np.std(lineZ.iloc[goodIdx])
    
    #assign line IDs
    lineIDs = physW.round(0)
    lineIDs = lineIDs.astype(int)
    
    #identify the line with highest SN
    strongestIdx = np.argmax(np.array(SNr))
    strongestID = lineIDs.iloc[strongestIdx]
                
    #if no lines are above threshold, use the line with highest SN
    if math.isnan(bestZ):
        strongestZ = lineZ.iloc[strongestIdx]
        bestZ = strongestZ
    
    objCols = ['objID','objFlag','splotZ','alfaZ', 'sigmaZ', 'nLines', 'bestID', 'redCoef',
             'redFlag', 'OIII/Hb', 'OII/Hb', 'NII/Ha', 'SII/Ha', 'NeIII/Hb', 'NeIII/OII']
    objDict = collections.OrderedDict(
                {'objID': [objID],
                'objFlag': [objFlag],
                'splotZ': format(z, '0.6f'),
                'alfaZ': format(bestZ, '0.6f'),
                'sigmaZ': format(bestZstd, '0.2e'),
                'nLines': nLines,
                'bestID': strongestID,
                'redCoef': np.nan,
                'redFlag': np.nan,
#                'redFlag': int(),
                'OIII/Hb': np.nan,
                'OII/Hb': np.nan,
                'NII/Ha': np.nan,
                'SII/Ha': np.nan,
                'NeIII/Hb': np.nan,
                'NeIII/OII': np.nan,
                })
    objDF = pd.DataFrame(objDict, columns=objCols)
    objDF.set_index(['objID'], inplace=True)
    objDF['objFlag'] = objDF['objFlag'].astype(int)
    objDF['nLines'] = objDF['nLines'].astype(int)
    objDF['bestID'] = objDF['bestID'].astype(int)
    
    lineCols = ['objID','lineID','observeW','physicalW','lineZ','flux','sigmaFlux','eqWidth','contAvg','fwhm']
    
    #god-awful, clunky re-formatting since there's just no good way to do this with pandas or numpy
    pd.options.mode.chained_assignment = None  # default='warn' - to mute the following warnings
    for i in range(nLines):
        observeW.iloc[i] = format( np.array(observeW)[i], '0.3f')
        lineZ.iloc[i] = format( np.array(lineZ)[i], '0.6f')
        flux.iloc[i] = format( np.array(flux)[i], '0.3e')
        sigmaFlux.iloc[i] = format( np.array(sigmaFlux)[i], '0.3e')
        eqWidth[i] = format( np.array(eqWidth)[i], '0.3f')
        contAvg[i] = format( np.array(contAvg)[i], '0.3e')
        fwhm.iloc[i] = format( np.array(fwhm)[i], '0.3f')
        
    lineDict = collections.OrderedDict(
                {'objID': [objID]*nLines,
                'lineID': lineIDs,
                'observeW': np.array(observeW),
                'physicalW': np.array(physW),
                'lineZ': np.array(lineZ),
                'flux': np.array(flux),
                'sigmaFlux': np.array(sigmaFlux),
                'eqWidth': eqWidth,
                'contAvg': contAvg,
                'fwhm': np.array(fwhm),
                })
                
    lineDF = pd.DataFrame(lineDict, columns=lineCols)
    lineDF.set_index(['objID'], inplace=True)
    return objDF, lineDF


def run1Dfrom2D(fpath, apNum):
    '''
    Run scopy on a spectrum specified by apNum, and process it with run1Dspec()
    '''
    pass


def callALFA(z, winSize, specFile, outPath):
    '''
    Internal function for run1Dspec. Calls the Automated-Line-Fitting-Algorithm
    (ALFA) from the command line using determined parameters.
    
    Args:
    z           redshift of the spectrum
    specFile    absolute path to the spectrum file
    outPath     absolute path to the output directory
    
    Returns:
    -----
    '''
    
    v = z * 299792	#km/s
    strongCat = ' optical_strongCustom.cat '
    deepCat = ' optical_deepCustom.cat '


    paramStr =  'alfa -ul  -n 0  -vg ' +str(v)+ ' -cw ' +str(winSize)+ \
                ' -o ' +outPath+ ' -sc' +strongCat+ '-dc' +deepCat+ ' '
    #-ul gets upper limit output; -n specifies normalization (none in this case)
    # -vg specifies velocity guess; -cw specifies the window size
    #-o specifies the path for output files; -sc/-dc specify line catalogs
    #See ALFA's manual for more info: >>> man alfa

    alfaCommand = paramStr + specFile
#    process = subprocess.call(alfaCommand, shell=True)
    process = subprocess.call(alfaCommand, shell=True, stdout=subprocess.PIPE)	#muted output


def promptLine(specFile, specType):
    '''
    Helper function - prompts user to measure a line in splot and to provide the
    identity of the measured line, finding both the rest and measured wavelengths.
    
    Provides shortcuts for certain lines, allowing for high precision for rest
    wavelengths.
    
    Args:
    specFile    absolute path to the spectrum file
    specType    type of spectrum to be measured; either 'sky' or 'science'
    
    Returns:
    restW       rest wavelength of one measured line in the spectrum
    userW       user-measured wavelength of one line in the spectrum, or:
    [empty]     only if user wants to skip line measuring (because e.g. the 
                spectrum has no good lines)
    '''
    ##Obtain measurement of a single line to get spectrum's redshift.
    #set up names and files
    specPath, _ = split(specFile)
    logFile = join(specPath, 'splot.log')
    #run splot to let user measure 1 line
    print('\n*****In the splot window, measure one line with a known rest ' \
            'wavelength, then press q. If there are no lines, just press q.')
    iraf.splot( images=specFile, save_file=logFile )

    ##Prompt user for the identity of the measured line
    while(True):    #continue until the user gives valid input
        if specType == 'science':   #user is measuring actual emission lines
            try:                    #ask for valid input
                key = eval(input(
                                '*****If you measured a line, enter:\n'
                                '0 for no lines - skip this spectrum\n'
                                '1 for H-alpha \n'
                                '2 for OIII \n'
                                '3 for OII \n'
                                'the wavelength of another line\n'
                                ))
                key = float(key)    #trigger a Name/ValueError if not a float
            except (NameError, ValueError, SyntaxError):
                print('Invalid input - please enter 0, 1, 2, 3, or a wavelength:\n')
                continue
            if key > 3 and key < 2000: continue
            #set the wavelength based on input key
            #it would be nice if Python had a switch-statement...
            if key == 0: return ('','',1)      #user wants to skip
            elif key == 1: restW = 6562.82      #H-alpha
            elif key == 2: restW = 5006.84      #OIII
            elif key == 3: restW = 3727.53      #OII blend
            else: restW = key
            break
        
        elif specType == 'sky':     #user is measuring sky emission lines
            try:                    #ask for valid input
                key = eval(input(
                                '*****If you measured a line, enter:\n'
                                '0 to skip\n'
                                '1 for Hg   4358\n'
                                '2 for Hg   5460\n'
                                '3 for O I  5577\n'
                                '4 for Na D 5891\n'
                                '5 for O I  6300\n'
                                '6 for O I  6363\n'
                                'the wavelength of another line\n'
                                ))
                key = float(key)    #trigger a Name/ValueError if not a float
            except (NameError, ValueError, SyntaxError):
                print('Invalid input - please enter a number key or a wavelength:\n')
                continue
            #set the wavelength based on input key
            if key == 0: return ('','',1)
            elif key == 1: restW = 4358.33      #Hg
            elif key == 2: restW = 5460.74      #Hg
            elif key == 3: restW = 5577.35      #OI
            elif key == 4: restW = 5891.94      #Na D
            elif key == 5: restW = 6300.23      #OI
            elif key == 6: restW = 6363.88      #OI
            else: restW = key
            break
    
    #read the measured wavelength from the splot log file
    entry = readSplotLog(logFile)
    userW = entry[0]
    skip = 0
    return (userW, restW, skip)


def classifySpectrum():
    
    validFlags = [0,1,2,3,4,5,6,17,18,19]
    print('\n*****Identify the spectrum with a flag:')
    while(True):     #loop until user gives valid input
        try: 
            objFlag = eval(input(
                            'To flag an object, enter:\n'
                            '0 for mystery\n'
                            '1 for Seyfert 1\n'
                            '2 for Seyfert 2\n'
                            '3 for star-forming\n'
                            '4 for LINER\n'
                            '5 for quasar\n'
                            '6 for HII Region\n'
                            '---\n'
                            '17 for false-detection (other)\n'
                            '18 for false-detection (galaxy)\n'
                            '19 for false-detection (star)\n'
                            ))
            objFlag = int(objFlag)    #trigger a Name/ValueError if not a float
        except (ValueError, NameError, SyntaxError):     #user gives non-integer input
            print('Invalid input - please enter a value from the list:\n')
            continue
        if objFlag not in validFlags: 
            print('Invalid input - please enter a value from the list:\n')
            continue
        break
    
    return objFlag


def readSplotLog(logFile):
    '''Helper function to read a splot log file, and extract the most recent entry'''
    #read splot.log to get the wavelength measured by user in splot
    cols = ['center','cont','flux','eqw','core','gfwhm','lfwhm']
    splotData = pd.read_table( logFile, names=cols, keep_default_na=False, sep='\s+' )
    #mark the invalid rows...
    splotData = splotData.apply(pd.to_numeric, errors='coerce')
    #...and get rid of them
    splotData = splotData.dropna(subset = ['center'])
    #get the most recent addition to the log file
    recentRow = splotData.iloc[-1][:]
    return recentRow

def sparseObjDF(objID, objFlag, splotZ=np.nan):
    '''
    User wants to skip a spectrum, so create a dataframe with the object ID and 
    object flag as the only global data.
    '''
    objCols = ['objID','objFlag','splotZ','alfaZ', 'sigmaZ', 'nLines', 'bestID', 'redCoef',
            'redFlag', 'OIII/Hb', 'OII/Hb', 'NII/Ha', 'SII/Ha', 'NeIII/Hb', 'NeIII/OII']
    objDict = collections.OrderedDict(
                {'objID': [objID],
                'objFlag': [objFlag],
                'splotZ': format(splotZ, '0.6f'),
                'alfaZ': np.nan,
                'sigmaZ': np.nan,
                'nLines': np.nan,
                'bestID': np.nan,
                'redCoef': np.nan,
                'redFlag': np.nan,
#                'redFlag': int(),
                'OIII/Hb': np.nan,
                'OII/Hb': np.nan,
                'NII/Ha': np.nan,
                'SII/Ha': np.nan,
                'NeIII/Hb': np.nan,
                'NeIII/OII': np.nan,
                })
    objDF = pd.DataFrame(objDict, columns=objCols)
    objDF.set_index(['objID'], inplace=True)
    return objDF

def sparseLineDF():
    '''
    User wants to skip a spectrum, so create a line dataframe with empty entries
    (for consistency, and to get proper output).
    '''
    lineCols = ['objID','lineID','observeW','physicalW','lineZ','flux','sigmaFlux','eqWidth','contAvg','fwhm']
    lineDF = pd.DataFrame(columns=lineCols)
    lineDF.set_index(['objID'], inplace=True)
    return lineDF

def runMultispec(fpath, skyFile=''):
    '''
    Identifies and measures the wavelengths, fluxes, and equivalent widths of 
    emission lines in a 2-dimensional set of spectra using Pyraf and ALFA.
    
    Prompts the user to measure and identify one line, via Pyraf's/Iraf's splot,
    in order to derive a redshift to be used in ALFA.
    
    Args:
    fpath       absolute path to a 2-dimensional sBased on your recent activity, I'd like to give you a special offer - a one month free upgrade to LinkedIn Premium.Based on your recent activity, I'd like to give you a special offer - a one month free upgrade to LinkedIn Premium.et of spectra .fits file
    skyFile     absolute path to a sky spectrum .fits files
    
    Returns:
    objDF           Pandas dataframe storing global data about each object. 
                    One entry per object.
    lineDF          Pandas dataframe storing data on emission lines in each 
                    spectrum. One entry per emission line for each spectrum.
    
    Output Files:
    /1dspectra      directory for 1D spectra generated from the input 2D spectrum
    1dspec.*.fits   a single 1D spectrum generated from the input 2D spectrum; 
    globalData.txt  csv output file containing global data about each object for 
                    each spectra that the user does not skip
    lineData.txt    same as above, but contains line data for each object
    '''

    datapath, basename = split(fpath)
    fname, __ = splitext(basename)
    
    #get the name of the data field, assuming a fname structure of "[NAME]_[EXTENSION]"
    fieldName = fname.split('_')[0]
    
    #directories for output files (create if nonexistent)
    specPath = join(datapath, '1dspectra')
    if not exists(specPath):
        makedirs(specPath)
    
    #prepare for sky line adjustments if desired
    if skyFile != '':
        try:
            print('\n*****Running a sky line correction. Displaying sky spectrum. ' \
                    'Follow the prompts and measure a sky line.')
            userW, restW, __ = promptLine(skyFile, 'sky')
            deltaW = userW-restW
        except: deltaW = 0
    else: deltaW = 0
    print('Sky line shift: ', round(deltaW,6), 'Angstroms' )
    
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
    
    #prompt user to say which spectrum to start on
    while(True):
        try:
            startSpec = eval(input('\nThere are ' +str(nSpec)+ ' spectra; press enter to start '
                        'at the first spectrum, or enter a number to start at a later spectrum:\n'))
            startSpec = int(startSpec)    #trigger a Name/ValueError if not an int
        except (NameError, ValueError):
            print('Invalid input - please enter an integer between ' +str(1)+ ' and ' +str(nSpec)+ ':\n')
            continue
        except SyntaxError:
            startSpec = 1
        if startSpec > nSpec or startSpec < 1:
            print('Invalid input - please enter an integer between ' +str(1)+ ' and ' +str(nSpec)+ ':\n')
            continue
        break

    
    ##Prepare to write to file. Each spectrum is processed and written with the 
    ##output file open, so processed data is still written even if the code 
    ##crashes/escapes. 
    outName1 = join(datapath, fieldName+'_globalData.txt')
    outName2 = join(datapath, fieldName+'_lineData.txt')
    
    
    objDF = pd.DataFrame()
    lineDF = pd.DataFrame()
    
    #ensure that data is still written despite errors
    with open(outName1, 'w+') as f1, open(outName2, 'w+') as f2:
        start = True
        for i in range(startSpec-1, nSpec):
#        for i in range(nSpec):
#        for i in range(39, 42):
        
            ##Run scopy to make a file for the current spectrum.

            objID = goodNames[i]
            #if different objects have the same objID, add 'a' to the latter duplicates
            #e.g. three entries of 2018 will become: '2018', '2018a', '2018aa'
            #This is redundantly called every loop, which could be cleaned up
            if (goodNames == objID).sum() > 1:
                nameIdx = np.argwhere(goodNames == objID)
                numDups = len(nameIdx)
                for x in range(numDups):
                    goodNames[nameIdx[x]] = objID + 'a'*x
            
            curAp = int(goodNums[i])
            pad_ap = format(curAp, '04d')
            specName = '1dspec.'+pad_ap+'.fits'
            specFile = join(specPath, specName)
            iraf.scopy(input=fpath, output=specFile, apertures=curAp, clobber='yes')
            
            #need to pause to let the file write
            sleep(1)
            
            ##Process the 1D spectrum.
            cur_objDF, cur_lineDF = run1Dspec(specFile, objID, deltaW)
        
            #make NaNs appear blank when viewing a file. But the blank spaces will still
            #be filled with NaNs when loaded back into Python (as desired)
            objDF = objDF.replace('nan','')
        
            objDF = objDF.append(cur_objDF)
            lineDF = lineDF.append(cur_lineDF)
            
            print('Writing line data for: ' +specFile+ '\n')
            if start:      #start of the file, so write headers
                f1.write( cur_objDF.to_csv(sep='\t', index=True) )
                f2.write( cur_lineDF.to_csv(sep='\t', index=True) )
                start = False
            else:
                f1.write( cur_objDF.to_csv(sep='\t', index=True, header=False) )
                f2.write( cur_lineDF.to_csv(sep='\t', index=True, header=False) )
    
    import completion; completion.thanksForPlaying()
    
    import lineRatios
    objDF, lineDF = lineRatios.redCorrRatios(objDF, lineDF)
    outName3 = join(datapath, fieldName+'_globalData_ratios.txt')
    objDF.to_csv(outName3, sep='\t', index=True)
    
    return objDF,lineDF
    
#    #scrape up all the files written by ALFA
#    allFiles = [f for f in listdir( join(datapath,savefolder) ) if f.endswith('.fits_fit') or f.endswith('.fits_lines')]



fpath = eval("input('Enter the full path of 2D spectrum fits file: ')")
skyFile = eval("input('If you want to apply sky line corrections, enter the full path to the sky spectrum. Otherwise, press enter:')")
#skyFile = '/home/bscousin/iraf/Team_SFACT/hadot055A/skyhadot055A_comb.fits'
#fpath = '/home/bscousin/iraf/Team_SFACT/hadot055A/hadot055A_comb_fin.ms.fits'

objDF, lineDF = runMultispec(fpath, skyFile)


