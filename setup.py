#!/usr/bin/python
from distutils.core import setup, Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

ext_modules = [Extension('chemeng.elementdata', ['src/chemeng/elementdata.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.thermodata', ['src/chemeng/thermodata.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.speciesdata', ['src/chemeng/speciesdata.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.antoinedata', ['src/chemeng/antoinedata.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.transportdata', ['src/chemeng/transportdata.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.components', ['src/chemeng/components.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.phase', ['src/chemeng/phase.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.reaction', ['src/chemeng/reaction.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.standarddefinitions', ['src/chemeng/standarddefinitions.pyx'], language='c++', include_dirs=['.']),
               Extension('chemeng.gibbminimizer', ['src/chemeng/gibbsminimizer.pyx'], language='c++', include_dirs=['.']),
               ]

setup(
    name="chemeng",
    version="0.1dev",
    author="M. Campbell Bannerman",
    author_email="m.bannerman@gmail.com",
    packages=['chemeng'],
    package_dir={'chemeng':'src/chemeng'},
    data_files=[#('PyChemEng', ['tests/chemeng.py']),
                ('PyChemEng/data', ['src/chemeng/data/antoine.inp',
                                    'src/chemeng/data/mass.mas03round.txt',
                                    'src/chemeng/data/isotopicCompositions.inp',
                                    'src/chemeng/data/thermo.inp'])],
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)
