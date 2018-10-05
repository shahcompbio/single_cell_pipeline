from setuptools import setup, find_packages
import versioneer


setup(
    name='single_cell',
    packages=find_packages(),
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Single cell pipeline',
    author='Andrew McPherson',
    author_email='andrew.mcpherson@gmail.com',
    entry_points={'console_scripts': ['single_cell = single_cell.run:main']},
    package_data={'':['scripts/*.py', 'scripts/*.R', 'scripts/*.npz', "config/*.yaml", "data/*"]}
)
