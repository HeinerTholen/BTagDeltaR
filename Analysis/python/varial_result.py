import itertools
import ROOT
import varial.analysis
import varial.tools


class FitResultCollector(varial.tools.Tool):
    def run(self):
        fit_chains = varial.analysis.lookup_children_names('../..')
        fitters = list(
            (name, varial.analysis.lookup(
                '../../%s/TemplateFitToolData' % name))
            for name in fit_chains
        )
        for name, res in fitters[:]:
            if res and 'TemplateFitTool' in res.name:
                res.name = name
            else:
                fitters.remove((name, res))
        self.result = sorted((res for _, res in fitters), key=lambda w: w.name)
        print self.result


class ScaleFactorsHisto(varial.tools.Tool):
    def __init__(self, name=None, with_cuts=False):
        super(ScaleFactorsHisto, self).__init__(name)
        self.with_cuts = with_cuts

    def run(self):
        inp = self.lookup('../FitResultCollector')
        inp = itertools.ifilter(lambda w:
                                (w.name[-4:] == 'Cuts') == self.with_cuts,
                                inp)


summary_chain = varial.tools.ToolChain(
    "Summary", [
        FitResultCollector(),
        ScaleFactorsHisto()
    ]
)