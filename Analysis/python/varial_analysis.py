import os

import ROOT
ROOT.gROOT.SetBatch()
import varial.main
import varial.tools

import ttdilep_samples

varial.tools.ToolChain(
    "ttdilep_analysis",
    [
        varial.tools.FwliteProxy(),
        varial.tools.FSStackPlotter()
    ]
)


varial.main.main(
    not_ask_execute=True,
    #fwlite_use_mp=False,
    samples=ttdilep_samples.smp_emu_mc,
    fwlite_executable=os.path.join(
        os.environ['CMSSW_BASE'],
        'src/BTagDeltaR/Analysis/python/vertexDR_worker.py',
    )
)