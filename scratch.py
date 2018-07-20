'''
Open and process WALF output??
'''

import pandas as pd
import numpy as np

from os.path import join
datapath = '/home/bscousin/iraf/Team_SFACT/hadot055A'
fname = 'hadot055A_comb_fin01.ms.txt'

import completion
completion.thanksForPlaying()


import numpy as np
def Round_To_n(x, n):
    return round(x, -int(np.floor(np.sign(x) * np.log10(abs(x)))) + n)
    


a = np.array([2.06333333e-17, 3.046666666E-17, 453.463, 5])
pd.set_option('display.float_format', '{:.3f}'.format)
df = pd.DataFrame(a)

df.to_csv('lines.csv', sep='\t', index=False)


objDF = pd.read_csv('globalData.txt', sep='\t', index_col=0)
lineDF = pd.read_csv('lineData.txt', sep='\t', index_col=0)

objDF.to_csv('globalData.txt',sep='\t', index=True) #good!!!
lineDF.to_csv('lineData.txt',sep='\t', index=True)
