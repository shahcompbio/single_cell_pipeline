'''
Created on 2013-09-28

@author: Andrew Roth
'''
import os
import random
import time

def get_ancestor_directory(path, level=1):
    '''
    Get the path of the directory a specified number of levels above the given path.
    
    >>> get_ancestor_directory('/foo/bar/some/where/my_file.txt', level=2)
    '/foo/bar/some'
    '''
    ancestor_dir = path
    
    for _ in range(level):
        ancestor_dir = os.path.dirname(ancestor_dir)
    
    return ancestor_dir

def make_directory(target_dir, mode=0775):
    '''
    Check if a directory exists and make it if not.
    
    For example, given /some/where make the folder /some/where. If /some does not exist, it will also be made.
    '''    
    i = 0
    
    try:
        old_umask = os.umask(0000)
        
        while not os.path.exists(target_dir):
            # Randomly sleep for a short random time so multiple simultaneous calls don't try to create the directory. 
            time.sleep(random.random() * 2)
            
            try:
                os.makedirs(target_dir, mode)
            
            except OSError:
                i += 1
                
                if i > 10:
                    raise
            
    finally:
        os.umask(old_umask)

def make_parent_directory(file_name, mode=0775):
    '''
    Given a file name, make the parent directory if it does not exist using make_directory.
    
    For example, given /some/where/foo.bar make the folder /some/where.
    '''
    parent_dir = os.path.dirname(file_name)
    
    make_directory(parent_dir, mode=mode)
    
if __name__ == '__main__':
    import doctest
     
    doctest.testmod()    