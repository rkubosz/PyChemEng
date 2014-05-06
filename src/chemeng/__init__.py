#!/usr/bin/env python
from chemeng.elementdata import elements
from chemeng.components import Components
from chemeng.speciesdata import speciesData,findSpeciesData
from chemeng.reaction import reaction
from chemeng.phase import IdealGasPhase, IncompressiblePhase
from chemeng.entropymaximiser import findEquilibrium

import chemeng.NASAdata
#import chemeng.chemkindata
