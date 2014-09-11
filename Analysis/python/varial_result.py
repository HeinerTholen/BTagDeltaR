import ctypes
import itertools
import pyparsing as p
import re
import os
import ROOT
import varial.analysis
import varial.dbio
import varial.tools


class TemplateFitPlots(varial.tools.Tool):
    def __init__(self, name=None, with_cuts=False):
        super(TemplateFitPlots, self).__init__(name)
        self.with_cuts = with_cuts

    def run(self):
        fit_chains = varial.analysis.lookup_children_names('../..')
        for name in fit_chains:
            if ('Cuts' in name) == self.with_cuts:
                os.system('cp %sVertexMass* %s' % (
                          varial.analysis.lookup_path(
                              '../../%s/TemplateFitToolData' % name),
                          self.cwd
                ))


class FitResultCollector(varial.tools.Tool):
    io = varial.dbio

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


class ScaleFactorsHisto(varial.tools.FSPlotter):
    pat = p.Suppress('FitChainSumfrom') \
          + p.Word(p.nums).setParseAction(lambda i: float(i[0])) \
          + p.Suppress('to') \
          + p.Word(p.nums).setParseAction(lambda i: float(i[0]))

    def __init__(self, name=None, with_cuts=False):
        super(ScaleFactorsHisto, self).__init__(name)
        self.with_cuts = with_cuts
        self.canvas_decorators = []

    def load_content(self):
        self.stream_content = itertools.ifilter(
            lambda w: (w.name[-4:] == 'Cuts') == self.with_cuts,
            self.lookup('../FitResultCollector')
        )

    def set_up_content(self):
        bins_cont = list(
            (self.pat.parseString(w.name)[:], (w.value, w.error, w.legend))
            for w in self.stream_content
        )
        n_bins = len(bins_cont)
        dbl_arr = ctypes.c_double * (n_bins + 1)
        bins = list(b[0]/20. for b, _ in bins_cont) + [bins_cont[-1][0][1]/20.]
        bins = dbl_arr(*bins)
        res = varial.wrappers.HistoWrapper(
            ROOT.TH1D(
                'ScaleFactors', ';#DeltaR(b, b); Data/MC',
                n_bins, bins
            ),
            draw_option='E1X1',
            analyzer=self.name,
        )
        for i in xrange(n_bins):
            val, err, leg = bins_cont[i][1]
            index = leg.index('B + B Vertex')
            res.histo.SetBinContent(i+1, val[index])
            res.histo.SetBinError(i+1, err[index])

        self.result = res
        self.stream_content = [[res]]


summary_chain = varial.tools.ToolChain(
    "Summary", [
        TemplateFitPlots(),
        TemplateFitPlots('TemplateFitPlotsCuts', True),
        FitResultCollector(),
        ScaleFactorsHisto(),
        ScaleFactorsHisto('ScaleFactorsHistoCuts', True),
    ]
)