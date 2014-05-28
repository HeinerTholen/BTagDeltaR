import os

import ROOT
ROOT.gROOT.SetBatch()
import varial.main
import varial.tools

import ttdilep_samples
samples = ttdilep_samples.smp_emu_mc


tc = varial.tools.ToolChain(
    "ttdilep_analysis",
    [
        varial.tools.FwliteProxy(),
        varial.tools.FSStackPlotter(),
        varial.tools.SimpleWebCreator(),
    ]
)


if __name__ == '__main__':
    varial.main.main(
        samples=samples,
        active_samples = list(s.name for s in samples),
        rootfile_postfixes=['.root', '.png'],
        fwlite_executable=os.path.join(
            os.environ['CMSSW_BASE'],
            'src/BTagDeltaR/Analysis/python/vertexDR_worker.py',
        ),
        toolchain=tc,
    )