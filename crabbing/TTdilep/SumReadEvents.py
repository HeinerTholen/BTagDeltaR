# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import glob
import itertools
import pyparsing as pp

# <codecell>

def gen_cat(f):
    with open(f) as fh:
        for line in fh.readlines():
            yield line
job_logs_strs = glob.glob("crab_0_140506_124239/res/CMSSW_*.stdout")
log_lines = itertools.chain.from_iterable(gen_cat(f) for f in job_logs_strs)

# <codecell>

parse_block = pp.Optional(pp.Suppress("TrigReport Events total =") 
                     + pp.Word(pp.nums).setParseAction(lambda t: int(t[0])) 
                     + pp.Suppress("passed ="))
parser = lambda s: parse_block.parseString(s).asList() or [0]

# <codecell>

sum_events = sum(parser(l)[0] for l in log_lines)
sum_events

# <codecell>


