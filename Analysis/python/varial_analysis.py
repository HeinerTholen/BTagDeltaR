import os

import ROOT
ROOT.gROOT.SetBatch()
import varial.main
import varial.tools
import varial.settings as s
import ttdilep_samples


s.web_target_dir = '/afs/desy.de/user/t/tholenhe/www/btagdr/ana/'
s.colors = ttdilep_samples.colors
s.rootfile_postfixes = ['.root', '.png']
s.fwlite_executable = os.path.join(
    os.environ['CMSSW_BASE'],
    'src/BTagDeltaR/Analysis/python/vertexDR_worker.py',
)


samples = ttdilep_samples.smp_emu_mc + ttdilep_samples.smp_emu_data
tc = varial.tools.ToolChain(
    "ttdilep_analysis",
    [
        varial.tools.FwliteProxy(),
        varial.tools.FSStackPlotter(keep_stacks_as_result=True),
        varial.tools.SimpleWebCreator(),
    ]
)


if __name__ == '__main__':
    varial.main.main(
        samples=samples,
        toolchain=tc,
        #max_num_processes=1,
    )
