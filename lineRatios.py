import pandas as pd
import numpy as np

def reddeningCoeff(objDF, lineDF):

    '''
    intrRatio   ratio between the intrinsic flux of the chosen line and 
                reference Hbeta line
    curveDiffBB    difference in reddening curve flux between the wavelengths of 
                the chosen line and reference Hbeta line
    '''
    pass

objDF = pd.read_csv('globalData.txt', sep='\t', index_col=0)
lineDF = pd.read_csv('lineData.txt', sep='\t', index_col=0)

Hbeta = 4861
Halpha = 6563
Hgamma = 4340

#load table of line IDs for later use
lineID_DF = pd.read_csv('lines.csv', sep='\t')
#reddening correction formula
reddeningCorr = lambda obsRatio, intrRatio, curveDiffB: -np.log10(obsRatio/intrRatio) / curveDiffB


objIDs = objDF.index


def findLineFlux(line):
    idx = np.where( objLines['lineID'].values == line)
    return float( objLines['flux'].iloc[idx] )

def findLineCurveDiff(line):
    ID_idx = np.where( lineID_DF['ID'].values == line)
    curLine = lineID_DF.iloc[ID_idx]
    return float(curLine['redCurve'])

def computeRatio(line1,line2,line3=None):
    line1Flux = findLineFlux(line1)
    if line3 != None:
        sumFlux = findLineFlux(line2) + findLineFlux(line3)
        line2Flux = sumFlux
    else:
        line2Flux = findLineFlux(line2)
    obsRatio = line2Flux / line1Flux
    
    curveDiffB1 = findLineCurveDiff(line1)
    curveDiffB2 = findLineCurveDiff(line2)
    
    corrRatio = obsRatio * 10**(redCoeff * (curveDiffB2-curveDiffB1) )
    return corrRatio

testIDs = ['1627','217','Dot']
for i in testIDs:
#for i in objIDs:
    #skip if the object has only 1 line (which triggers errors if I don't skip it)
    if objDF.loc[i]['nLines'] == 1: continue
    
    #look at only the current object's lines
    objLines = lineDF.loc[i]
    #Halpha detected
    if Halpha in objLines['lineID'].values and Hbeta in objLines['lineID'].values:
        line = Halpha
        obsFlux = findLineFlux(line)
        intrRatio = 2.86        #Halpha:Hbeta
        flag = 1
            
    #OIII detected
    elif Hgamma in objLines['lineID'].values and Hbeta in objLines['lineID'].values:
        line = Hgamma
        obsFlux = findLineFlux(line)
        intrRatio = 0.47        #Hgamma:Hbeta
        flag = 2
    else: continue #cry me a river

    #look up values in line ID table
    curveDiffB = findLineCurveDiff(line)

    HbetaFlux = findLineFlux(Hbeta)
    obsRatio = float(obsFlux / HbetaFlux)

    redCoeff = reddeningCorr(obsRatio, intrRatio, curveDiffB)

    objDF.loc[i, 'redCoeff'] = redCoeff
#    objDF.loc[i, 'redFlag'] = flag
    
    ##Now compute line ratios
    Hbeta = 4861
    Halpha = 6563
    Hgamma = 4340
    NII = 6584
    SII_1 = 6717
    SII_2 = 6731
    OII_1 = 3726
    OII_2 = 3729
    OIII = 5007
    
    NeIII = 3869
    if Halpha in objLines['lineID'].values:
        line1 = Halpha
        if NII in objLines['lineID'].values:
            line2 = NII
            corrRatio = computeRatio(line1,line2)
            objDF.loc[i, 'NII/Ha'] = corrRatio
#            import pdb; pdb.set_trace()
        if SII_1 in objLines['lineID'].values and SII_2 in objLines['lineID'].values:
            line2 = SII_1
            line3 = SII_2
            corrRatio = computeRatio(line1,line2,line3)
            objDF.loc[i, 'SII/Ha'] = corrRatio
            
    if Hbeta in objLines['lineID'].values:
        line1 = Hbeta
        if OIII in objLines['lineID'].values:
            line2 = OIII
            corrRatio = computeRatio(line1,line2)
            objDF.loc[i, 'OIII/Hb'] = corrRatio
            
        if OII_1 in objLines['lineID'].values and OII_2 in objLines['lineID'].values:
            line2 = OII_1
            line3 = OII_2
            corrRatio = computeRatio(line1,line2,line3)
            objDF.loc[i, 'OII/Hb'] = corrRatio
                        
        if NeIII in objLines['lineID'].values:
            line2 = NeIII
            corrRatio = computeRatio(line1,line2)
            objDF.loc[i, 'NeIII/Hb'] = corrRatio
    
    if NeIII in objLines['lineID'].values and OII_1 in objLines['lineID'].values and OII_2 in objLines['lineID'].values:
        #line1 = OII
        #line2 = NeIII
        line1 = NeIII
        line2 = OII_1
        line3 = OII_2
        corrRatio = 1. /computeRatio(line1,line2)
        objDF.loc[i, 'NeIII/OII'] = corrRatio
    import pdb; pdb.set_trace()

        
def lineRatios(objDF, lineDF):
    pass
    '''
    'OIII/Hb', 'OII/Hb', 'NII/Ha', 'SII/Ha', 'NeIII/Hb', 'NeIII/OII'
    '''


'''
    if Halpha in objLines['lineID'].values:
        HalphaFlux = findLineFlux(Halpha)
        curveDiffB1 = findLineCurveDiff(Halpha)
        if NII in objLines['lineID'].values:
            line = NII
            obsFlux = findLineFlux(line)
            curveDiffB2 = findLineCurveDiff(line)
            obsRatio = obsFlux/HalphaFlux
            corrRatio = obsRatio * 10**(redCoeff * curveDiffB2-curveDiffB1)
            objDF.loc[i, 'NII/Ha'] = corrRatio
            
        if SII in objLines['lineID'].values:
            line = SII
            obsFlux = findLineFlux(line)
            objDF.loc[i, 'SII/Ha'] = obsFlux/HalphaFlux
            
    if Hbeta in objLines['lineID'].values:
        HbetaFlux = findLineFlux(Hbeta)
        if OIII in objLines['lineID'].values:
            line = OIII
            obsFlux = findLineFlux(line)
            objDF.loc[i, 'OIII/Hb'] = obsFlux/HbetaFlux
            
        if OII in objLines['lineID'].values:
            line = OII
            obsFlux = findLineFlux(line)
            objDF.loc[i, 'OII/Hb'] = obsFlux/HbetaFlux
            
        if NeIII in objLines['lineID'].values:
            line = NeIII
            obsFlux = findLineFlux(line)
            objDF.loc[i, 'NeIII/Hb'] = obsFlux/HbetaFlux
    
    if NeIII in objLines['lineID'].values and OII in objLines['lineID'].values:
        line = NeIII
        obsFlux1 = findLineFlux(line)
        line = OII
        obsFlux2 = findLineFlux(line)
        objDF.loc[i, 'NeIII/OII'] = obsFlux1/obsFlux2
#    import pdb; pdb.set_trace()

'''




'''
Halpha = 6563
Hgamma = 4340
Hbeta = 4861
def huh(objDF, lineDF, wave1, wave2, waveRef):
    if wave1 in lineDF['lineID'].values and waveRef in lineDF['lineID'].values:
        #Halpha
        intrRatio = 2.86
        idx = np.where( lineDF['lineID'].values == wave1)
    
        
    elif wave2 in lineDF['lineID'].values and waveRef in lineDF['lineID'].values:
        #Hgamma
        intrRatio = 0.47
        idx = np.where( lineDF['lineID'].values == wave2)
        
    obsFlux = lineDF['flux'].iloc[idx]
    
    ref_idx = np.where( lineDF['lineID'].values == waveRef)
    refFlux = lineDF['flux'].iloc[ref_idx]
    

    
    return curveDiffB
'''
