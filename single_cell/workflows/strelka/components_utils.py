'''
Created on Nov 21, 2015

@author: Andrew Roth
'''
import os
import random
import time
import errno


def find(name, path):
    for root, _, files in os.walk(path):
        if name in files:
            return os.path.join(root, name)


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


def flatten_input(files):
    if type(files) == dict:
        parsed_files = [files[x] for x in sorted(files)]
    elif type(files) == str:
        parsed_files = [files, ]
    else:
        parsed_files = []
        for x in files:
            if type(x) == dict:
                parsed_files.extend([x[y] for y in sorted(x)])
            else:
                parsed_files.append(x)
    return parsed_files


def remove(filename):
    '''
    Remove a file that may not exist
    '''
    try:
        os.remove(filename)
    except OSError as e:
        if e.errno != errno.ENOENT:
            raise


def symlink(filename, link_name=None, link_directory=None):
    '''
    Create a symlink, with additional options for flexibility,

    Args:
        filename (str): file to link to

    KwArgs:
        link_name (str): base name of the link, defaults to same as link to
        link_directory (str): directory of the, defaults to directory of link to

    '''
    if link_name is None:
        link_name = os.path.basename(filename)
    if link_directory is None:
        link_directory = os.getcwd()
    link_filename = os.path.join(link_directory, link_name)
    remove(link_filename)
    filename = os.path.abspath(filename)
    os.symlink(filename, link_filename)
    return link_filename

if __name__ == '__main__':
    import doctest

    doctest.testmod()
