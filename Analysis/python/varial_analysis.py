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
active_samples = list(s.name for s in samples)
active_samples.remove('TTbar')

tc = varial.tools.ToolChain(
    "ttdilep_analysis",
    [
        varial.tools.FwliteProxy(fwlite_exe),
        varial.tools.ZipTool('ttdilep_analysis/FwliteProxy'),
        varial.tools.CopyTool(os.path.join(os.environ['HOME'], 'www/btagdr/ana/'), name="ZipFileCopyTool",),
        varial_plotters.jet_plots,
        varial_plotters.SampleNormalizer(),
    #] + varial_plotters.chains + [
        varial_fitter.fitter_plots,
        varial_fitter.fitter_chain,
        varial_fitter.fitter_chain_sum,
        varial.tools.SimpleWebCreator(),
        varial.tools.CopyTool(
            os.path.join(os.environ['HOME'], 'www/btagdr/ana/')),
    ]
)


def main():
    varial.main.main(
        samples=samples,
        active_samples=active_samples,
        toolchain=tc,
    )


if __name__ == '__main__':
    main()

