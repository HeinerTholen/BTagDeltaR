import os

import ROOT
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
import varial.main
import varial.tools
import varial.settings as s

import ttdilep_samples
import varial_plotters
import varial_fitter

s.rootfile_postfixes = ['.root', '.png']
fwlite_exe = os.path.join(
    os.environ['CMSSW_BASE'],
    'src/BTagDeltaR/Analysis/python/worker_vertexDR.py',
)

samples = ttdilep_samples.smp_emu_mc + ttdilep_samples.smp_emu_data
tc = varial.tools.ToolChain(
    "ttdilep_analysis",
    [
        varial.tools.FwliteProxy(fwlite_exe),
        varial.tools.ZipTool('ttdilep_analysis/FwliteProxy'),
        varial_plotters.DaNormalizer(),
        varial_plotters.chain_ivf_merged_filt,
        varial_plotters.chain_ivf_b2c_merged,
        varial_plotters.chain_ivf_b2c_merged_filt,
        varial_fitter.fitter_chain,
        varial.tools.SimpleWebCreator(),
        varial.tools.CopyTool(
            os.path.join(os.environ['HOME'], 'www/btagdr/ana/')),
    ]
)


def main():
    varial.main.main(
        samples=samples,
        toolchain=tc,
    )


if __name__ == '__main__':
    main()

