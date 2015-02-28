#!/usr/bin/python
from distutils.core import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

ext_modules = [Extension('chemeng.elementdata', ['src/chemeng/elementdata.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.NASAdata', ['src/chemeng/NASAdata.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.cementdata', ['src/chemeng/cementdata.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.chemkindata', ['src/chemeng/chemkindata.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.speciesdata', ['src/chemeng/speciesdata.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.antoinedata', ['src/chemeng/antoinedata.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.transportdata', ['src/chemeng/transportdata.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.components', ['src/chemeng/components.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.ponchonsavarit', ['src/chemeng/PonchonSavarit.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.phase', ['src/chemeng/phase.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.standarddefinitions', ['src/chemeng/standarddefinitions.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.entropymaximiser', ['src/chemeng/entropymaximiser.pyx'], language='c++', include_dirs=['.']),
               ]

setup(
    name="chemeng",
    version="0.1dev",
    author="M. Campbell Bannerman",
    author_email="m.bannerman@gmail.com",
    packages=['chemeng'],
    package_dir={'chemeng':'src/chemeng'},
    package_data={'chemeng' : ['data/antoine.inp',
                               'data/mass.mas03round.txt',
                               'data/isotopicCompositions.inp',
                               'data/NASA_CEA.inp',
                               'data/NEWNASA.TXT',
                               'data/Cement.csv',
                               'data/Cement2.csv',
                               'data/BurcatCHEMKIN.DAT',
                               'data/NistData.csv',
                               'data/Cement_New_Tests.csv',
                               'data/Cement_Therm_New2.csv'
                           ]},
    cmdclass = {'build_ext': build_ext},
    py_modules = ['chemeng.config'],
    ext_modules = ext_modules
)
