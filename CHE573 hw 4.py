# -*- coding: utf-8 -*-
"""

Created on Sun Oct  3 16:00:56 2021

@author: Travis J Czechorski

Github: https://github.com/tjczec01

Email:   travis.czechorski@maine.edu
         tjczec01@gmail.com
        
Advisor: thomas.schwartz@maine.edu
       
Github:  https://github.com/tjczec01 
         https://github.com/TravisCzechorskiUMaine

"""

import pandas as pd
import math as mt

df = pd.read_excel(r'excilefilepath.xlsx')
energy_list = df['-205.03321714472892'].tolist()
# print(energy_list)

def std_dev(pop_list):
    N = len(pop_list)
    AVG = sum(pop_list)/len(pop_list)
    subtracted_values = [abs(i - AVG) for i in pop_list]
    numerator = sum(subtracted_values)
    final_value = mt.sqrt(numerator/N)
    return final_value

answer = std_dev(energy_list)
print(answer)
    
