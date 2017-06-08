from distutils.core import setup

setup(
      name='Pipelines',
      version='0.2.8',
      description='Useful tools for working with Ruffus pipelines.',
      author='Andrew Roth',
      author_email='andrewjlroth@gmail.com',   
      packages=[ 
                'pipelines'
                ],
      data_files = [('', ['LICENSE.txt', 'README.md'])]
     )
