#!/usr/bin/env python
import os,sys,pdb,scipy,glob
from pylab import *
from strolger_util import util as u
from strolger_util import rates_z as rz



p = (-511.1153201481045, 53.577081653663541, -2.0012061989355914)

times = arange(0, 13.6, 0.05)

tmp = rz.dtdfunc(times, *p)

pdb.set_trace()
