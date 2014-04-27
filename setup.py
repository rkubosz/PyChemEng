#!/usr/bin/python
from distutils.core import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

ext_modules = [Extension('chemeng.elementdata', ['chemeng/elementdata.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.thermodata', ['chemeng/thermodata.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.speciesdata', ['chemeng/speciesdata.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.antoinedata', ['chemeng/antoinedata.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.transportdata', ['chemeng/transportdata.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.components', ['chemeng/components.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.phase', ['chemeng/phase.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.reaction', ['chemeng/reaction.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.standarddefinitions', ['chemeng/standarddefinitions.pyx'], language='c++', include_dirs=['.']),
               ]

setup(
    name="chemeng",
    version="0.1dev",
    author="M. Campbell Bannerman",
    author_email="m.bannerman@gmail.com",
    packages=['chemeng'],
    data_files=[#('PyChemEng', ['tests/chemeng.py']),
                ('PyChemEng/data', ['chemeng/data/antoine.inp',
                                    'chemeng/data/mass.mas03round.txt',
                                    'chemeng/data/isotopicCompositions.inp',
                                    'chemeng/data/thermo.inp'])],
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
