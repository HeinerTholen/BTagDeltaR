import os

import ROOT
ROOT.gROOT.SetBatch()
ROOT.gROOT.ProcessLine('gErrorIgnoreLevel = kError;')
import varial.main
import varial.tools
import varial.settings as s

import ttdilep_samples
import varial_plotters

s.web_target_dir = '/afs/desy.de/user/t/tholenhe/www/btagdr/ana/'
s.rootfile_postfixes = ['.root', '.png']
s.fwlite_executable = os.path.join(
    os.environ['CMSSW_BASE'],
    'src/BTagDeltaR/Analysis/python/worker_vertexDR.py',
)

samples = ttdilep_samples.smp_emu_mc + ttdilep_samples.smp_emu_data
tc = varial.tools.ToolChain(
    "ttdilep_analysis",
    [
        varial.tools.FwliteProxy(),
        varial_plotters.DaNormalizer(),
#        varial_plotters.chain_ivf_merged,
        varial_plotters.chain_ivf_merged_filt,
        varial_plotters.chain_ivf_b2c_merged,
        varial_plotters.chain_ivf_b2c_merged_filt,
#        varial_plotters.chain_ivf_b2c_merged_filt_cov,
        varial.tools.SimpleWebCreator(),
    ]
)


def main():
    varial.main.main(
        samples=samples,
        toolchain=tc,
    )


if __name__ == '__main__':
    main()

