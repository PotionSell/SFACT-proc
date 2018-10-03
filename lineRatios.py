import math
import pandas as pd
import numpy as np

nullEntry = -1
numDec = 4      #number of decimal places, for formatting

def redCorrRatios(objDF, lineDF):

    '''
    To be updated...
    intrRatio   ratio between the intrinsic flux of the chosen line and
                reference Hbeta line
    curveDiffBB    difference in reddening curve flux between the wavelengths of
                the chosen line and reference Hbeta line

    If the computed reddening coefficient is negative, then it is recorded in the
    output but is NOT used for correcting the line ratios.
    '''
    nullVal = -9.999999

#    objDF = pd.read_csv('globalData00.txt', sep='\t', index_col=0)
#    lineDF = pd.read_csv('lineData00.txt', sep='\t', index_col=0)

    ##Internal functions:
    #reddening correction formula
    reddeningCorr = lambda obsRatio, intrRatio, curveDiffB: -np.log10(obsRatio/intrRatio) / curveDiffB

    def findLineFlux(line):
        idx = np.where( objLines['lineID'].values == line)
        return float( objLines['flux'].iloc[idx] )

    def findLineCurveDiff(line):
        ID_idx = np.where( lineID_DF['ID'].values == line)
        curLine = lineID_DF.iloc[ID_idx]
        return float(curLine['redCurve'])

    def computeRatio(flag,line1,line2,line3=None):
        line1Flux = findLineFlux(line1)
        if line3 != None:
            sumFlux = findLineFlux(line2) + findLineFlux(line3)
            line2Flux = sumFlux
        else:
            line2Flux = findLineFlux(line2)

        #multiply by a flat value, so that machine precision uncertainties for
        #small numbers is avoided
        scale = 10E6
        line2Flux = line2Flux*scale
        line1Flux = line1Flux*scale
        #the scale divides out in the ratio
        obsRatio = line2Flux / line1Flux
#        import pdb; pdb.set_trace()
        if flag == nullEntry:
            corrRatio = obsRatio
        else:
            curveDiffB1 = findLineCurveDiff(line1)
            curveDiffB2 = findLineCurveDiff(line2)
            #if there's a negative redCoef, then record it but do NOT use it for corrections
            if redCoef < 0:
                return obsRatio
            corrRatio = obsRatio * 10**(redCoef * (curveDiffB2-curveDiffB1) )
        return corrRatio


    Hbeta = 4861
    Halpha = 6563
    Hgamma = 4340

    #load table of line IDs for later use
    lineID_DF = pd.read_csv('lines.csv', sep='\t')

    objIDs = objDF.index

    #testIDs = ['1627','217','Dot']
    #for i in testIDs:

    for i in objIDs:
        #look at only the current object's lines
        objLines = lineDF.loc[[i]]      #extra [ ] creates a 1D df, instead of a series
        nLines = len(objLines)
        #skip if the object has zero or one line (which triggers errors if I don't skip it)
        if nLines == 1 or np.isnan(nLines):
            flag = nullEntry
            redCoef = nullEntry
            objDF.loc[i, 'redCoef'] = redCoef
            objDF.loc[i, 'redFlag'] = flag
            continue

        #Halpha detected
        if Halpha in objLines['lineID'].values and Hbeta in objLines['lineID'].values:
            line = Halpha
            obsFlux = findLineFlux(line)
            intrRatio = 2.86        #Halpha:Hbeta
            flag = 1

            curveDiffB = findLineCurveDiff(line)

            HbetaFlux = findLineFlux(Hbeta)
            obsRatio = float(obsFlux / HbetaFlux)

            redCoef = reddeningCorr(obsRatio, intrRatio, curveDiffB)

        #OIII detected
        elif Hgamma in objLines['lineID'].values and Hbeta in objLines['lineID'].values:
            line = Hgamma
            obsFlux = findLineFlux(line)
            intrRatio = 0.47        #Hgamma:Hbeta
            flag = 2

            curveDiffB = findLineCurveDiff(line)

            HbetaFlux = findLineFlux(Hbeta)
            obsRatio = float(obsFlux / HbetaFlux)

            redCoef = reddeningCorr(obsRatio, intrRatio, curveDiffB)
        else:   #cry me a river
            flag = nullEntry
            redCoef = nullEntry

        #if a coefficient has already been stored, read it in instead
        if not math.isnan(objDF.loc[i, 'redCoef']):
            redCoef = objDF.loc[i, 'redCoef']
            flag = objDF.loc[i, 'redFlag']
        else:
            objDF.loc[i, 'redCoef'] = np.round( redCoef, numDec )
            objDF.loc[i, 'redFlag'] = np.round( flag, numDec )

        ##Compute line ratios
        NII = 6584
        SII_1 = 6716
        SII_2 = 6731
        OII_1 = 3726
        OII_2 = 3729
        OIII = 5007
        NeIII = 3869



        if Halpha in objLines['lineID'].values:
            line1 = Halpha
            if NII in objLines['lineID'].values:
                line2 = NII
                corrRatio = computeRatio(flag,line1,line2)
                objDF.loc[i, 'NII/Ha'] = np.round(np.log10(corrRatio), numDec)
            if SII_1 in objLines['lineID'].values and SII_2 in objLines['lineID'].values:
                line2 = SII_1
                line3 = SII_2
                corrRatio = computeRatio(flag,line1,line2,line3)
                objDF.loc[i, 'SII/Ha'] = np.round(np.log10(corrRatio), numDec)

        if Hbeta in objLines['lineID'].values:
            line1 = Hbeta
            if OIII in objLines['lineID'].values:
                line2 = OIII
                corrRatio = computeRatio(flag,line1,line2)
                objDF.loc[i, 'OIII/Hb'] = np.round(np.log10(corrRatio), numDec)
            if OII_1 in objLines['lineID'].values and OII_2 in objLines['lineID'].values:
                line2 = OII_1
                line3 = OII_2
                corrRatio = computeRatio(flag,line1,line2,line3)
                objDF.loc[i, 'OII/Hb'] = np.round(np.log10(corrRatio), numDec)
            if NeIII in objLines['lineID'].values:
                line2 = NeIII
                corrRatio = computeRatio(flag,line1,line2)
                objDF.loc[i, 'NeIII/Hb'] = np.round(np.log10(corrRatio), numDec)

        if NeIII in objLines['lineID'].values and OII_1 in objLines['lineID'].values and OII_2 in objLines['lineID'].values:
            #line1 = OII
            #line2 = NeIII
            line1 = NeIII
            line2 = OII_1
            line3 = OII_2
            corrRatio = 1. / computeRatio(flag,line1,line2,line3)
            objDF.loc[i, 'NeIII/OII'] = np.round(np.log10(corrRatio), numDec)
    objDF['redFlag'] = objDF['redFlag'].astype(int)

    #make NaNs appear blank when viewing a file. But the blank spaces will still
    #be filled with NaNs when loaded back into Python (as desired)
#    objDF = objDF.replace('nan','')
    return objDF, lineDF

#objDF.to_csv('globalData.txt',sep='\t', index=True)
#lineDF.to_csv('lineData.txt',sep='\t', index=True)
