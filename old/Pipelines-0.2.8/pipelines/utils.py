'''
Created on 2014-04-05

@author: Andrew Roth
'''
from collections import OrderedDict

import os

def get_sample_information(file_name):
    '''
    Parse a file name into a dictionary. Key-value pairs are delimited by a ".", and the key is the value before the
    first '_'.
    
    Examples: 
    - /foo/bar/key_value.tsv
    - /foo/bar/long-key_long-value.key2_value2.xml
    - /foo/bar/key_value_value.xls
    '''
    file_name = os.path.basename(file_name)
    
    information = OrderedDict()
    
    for x in file_name.split('.'):
        parts = x.split('_')
        
        if len(parts) > 1:
            information[parts[0]] = '_'.join(parts[1:])
    
    return information