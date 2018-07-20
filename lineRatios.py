def reddeningCoeff(objDF, lineDF):
    pass

    '''
    intrRatio   ratio between the intrinsic flux of the chosen line and 
                reference Hbeta line
    curveDiff    difference in reddening curve flux between the wavelengths of 
                the chosen line and reference Hbeta line
    '''

objDF = pd.read_csv('globalData.txt', sep='\t', index_col=0)
lineDF = pd.read_csv('lineData.txt', sep='\t', index_col=0)

Hbeta = 4861
Halpha = 6563
Hgamma = 4340

#load table of line IDs for later use
lineID_DF = pd.read_csv('lines.csv', sep='\t')
#reddening correction formula
reddeningCorr = lambda obsRatio, intrRatio, curveDiff: -np.log10(obsRatio/intrRatio) / curveDiff


objIDs = objDF.index


def findLineFlux(objLines, line):
    idx = np.where( objLines['lineID'].values == line)
    return float( objLines['flux'].iloc[idx] )


for i in objIDs:
    #look at only the current object's lines
    objLines = lineDF.loc[i]
    #Halpha detected
    if Halpha in objLines['lineID'].values and Hbeta in objLines['lineID'].values:
        line = Halpha
        obsFlux = findLineFlux(objLines, line)
        intrRatio = 2.86        #Halpha:Hbeta
            
    #OIII detected
    elif Hgamma in objLines['lineID'].values and Hbeta in objLines['lineID'].values:
        line = Hgamma
        obsFlux = findLineFlux(objLines, line)
        intrRatio = 0.47        #Hgamma:Hbeta
    else: continue #cry me a river

    #look up values in line ID table
    ID_idx = np.where( lineID_DF['ID'].values == line)
    curLine = lineID_DF.iloc[ID_idx]
    curveDiff = float(curLine['redCurve'])

    HbetaFlux = findLineFlux(objLines, Hbeta)
    obsRatio = float(obsFlux / HbetaFlux)

    redCoeff = reddeningCorr(obsRatio, intrRatio, curveDiff)

    objDF.loc[i, 'redCoeff'] = redCoeff

    
    ##Now compute line ratios
    Hbeta = 4861
    Halpha = 6563
    Hgamma = 4340
    NII = 6548
    SII = 4072
    OIII = 5007
    OII = 3727      #this is blended - need  to account for the rough flux!!
    NeIII = 3970
    if Halpha in objLines['lineID'].values:
        HalphaFlux = findLineFlux(objLines, Halpha)
        if NII in objLines['lineID'].values:
            line = NII
            obsFlux = findLineFlux(objLines, line)
            objDF.loc[i, 'NII/Ha'] = obsFlux/HalphaFlux
            
        if SII in objLines['lineID'].values:
            line = SII
            obsFlux = findLineFlux(objLines, line)
            objDF.loc[i, 'SII/Ha'] = obsFlux/HalphaFlux
            
    if Hbeta in objLines['lineID'].values:
        HbetaFlux = findLineFlux(objLines, Hbeta)
        if OIII in objLines['lineID'].values:
            line = OIII
            obsFlux = findLineFlux(objLines, line)
            objDF.loc[i, 'OIII/Hb'] = obsFlux/HbetaFlux
            
        if OII in objLines['lineID'].values:
            line = OII
            obsFlux = findLineFlux(objLines, line)
            objDF.loc[i, 'OII/Hb'] = obsFlux/HbetaFlux
            
        if NeIII in objLines['lineID'].values:
            line = NeIII
            obsFlux = findLineFlux(objLines, line)
            objDF.loc[i, 'NeIII/Hb'] = obsFlux/HbetaFlux
    
    if NeIII in objLines['lineID'].values and OII in objLines['lineID'].values:
        line = NeIII
        obsFlux1 = findLineFlux(objLines, line)
        line = OII
        obsFlux2 = findLineFlux(objLines, line)
        objDF.loc[i, 'NeIII/OII'] = obsFlux1/obsFlux2
#    import pdb; pdb.set_trace()

        
def lineRatios(objDF, lineDF):
    pass
    '''
    'OIII/Hb', 'OII/Hb', 'NII/Ha', 'SII/Ha', 'NeIII/Hb', 'NeIII/OII'
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
    

    
    return curveDiff
'''
