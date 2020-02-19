#!/usr/bin/env python
from __future__ import division
import sys
import os
import argparse

# Path hack.
import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.insert(0,currentdir + "/parameter_definitions")


from chosen_parameter_sets import getParameterSets

for paramSet in getParameterSets():
    paramSet.analyze()
    print "================================================================================"
    print "%.80s" % ("===== Round5 %s ============================================================================" % paramSet.name)
    print "================================================================================"
    print
    print paramSet
    print
