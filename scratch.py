'''
Open and process WALF output??

huh??
'''

#    #scrape up all the files written by ALFA
#    allFiles = [f for f in listdir( join(datapath,savefolder) ) if f.endswith('.fits_fit') or f.endswith('.fits_lines')]

import pandas as pd
import numpy as np

from os.path import join
datapath = '/home/bscousin/iraf/Team_SFACT/hadot055A'
fname = 'hadot055A_comb_fin01.ms.txt'

import completion
completion.thanksForPlaying()


###formatting shenanigans
#trying a function
import numpy as np
def Round_To_n(x, n):
    return round(x, -int(np.floor(np.sign(x) * np.log10(abs(x)))) + n)
#

#tests
a = np.array([2.06333333e-17, 3.046666666E-17, 453.463, 5])
pd.set_option('display.float_format', '{:.3f}'.format)
df = pd.DataFrame(a)

df.to_csv('lines.csv', sep='\t', index=False)


objID = objID.ljust(4)  #pad for output formatting purposes
floatFormatter= lambda x: format(x, '6.3E')
np.set_printoptions(formatter={'float': floatFormatter})

###

###i/o with dataframes
datapath = '/home/bscousin/iraf/Team_SFACT/hadot055A'
objPath = join(datapath, 'globalDataJJS.txt')
objDF = pd.read_csv(objPath, sep='\t', index_col=0)
linePath = join(datapath, 'lineDataJJS.txt')
lineDF = pd.read_csv(linePath, sep='\t', index_col=0)

objDF.to_csv('globalData.txt',sep='\t', index=True) #good!!!
lineDF.to_csv('lineData.txt',sep='\t', index=True)
###

import pickle
pickle.dump(fig, open(figname, 'wb'))
pickle.load(open(figname, 'rb'))
