from setuptools import setup, find_packages, Extension

setup(
    name='single_cell_nextseq',
    packages=find_packages(),
    description='Single cell pipeline for nextseq',
    author='Andrew McPherson',
    author_email='andrew.mcpherson@gmail.com',
    entry_points={'console_scripts': ['single_cell_nextseq = single_cell_nextseq.run:main']},
    package_data={'':['scripts/*']}
)
