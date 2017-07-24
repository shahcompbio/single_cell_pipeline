from setuptools import setup, find_packages, Extension

setup(
    name='single_cell',
    packages=find_packages(),
    description='Single cell pipeline',
    author='Andrew McPherson',
    author_email='andrew.mcpherson@gmail.com',
    entry_points={'console_scripts': ['single_cell = single_cell.run:main']},
    package_data={'':['scripts/*']}
)
