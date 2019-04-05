#!/usr/bin/env python
from __future__ import division
import sys
import os
import argparse
from parameter_sets import getParameterSets

for paramSet in getParameterSets():
    paramSet.analyze()
    print "================================================================================"
    print "%.80s" % ("===== Round5 %s ============================================================================" % paramSet.name)
    print "================================================================================"
    print
    print paramSet
    print
