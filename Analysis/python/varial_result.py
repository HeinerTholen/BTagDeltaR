import ctypes
import itertools
import pyparsing as p
import re
import os
import ROOT
import varial.analysis
import varial.dbio
import varial.tools
import varial.util


class TemplateFitPlots(varial.tools.Tool):
    """Copy template fit plots into one directory."""
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


class TemplatePlots(TemplateFitPlots):
    def run(self):
        fit_chains = varial.analysis.lookup_children_names('../..')
        if self.with_cuts:
            cp_str = 'cp %sIvfB2cMergedCuts_VertexMass* %s'
        else:
            cp_str = 'cp %sIvfB2cMerged_VertexMass* %s'
        for name in fit_chains:
            if ('Cuts' in name) == self.with_cuts:
                sys_call = cp_str % (
                          varial.analysis.lookup_path(
                              '../../%s/TemplatesAndFittedHisto' % name),
                          self.cwd
                )
                self.message("INFO " + sys_call)
                os.system(sys_call)


class FitResultCollector(varial.tools.Tool):
    """Collect fit results (numbers only)."""
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


######################################################### ScaleFactorsHisto ###
class IvfEfficiencyLine(varial.util.Decorator):
    def do_final_cosmetics(self):
        self.decoratee.do_final_cosmetics()
        y_val = .92**2
        line = ROOT.TLine()
        line.SetLineWidth(2)
        line.SetLineColor(ROOT.kSpring - 4)
        line.DrawLine(0., y_val, 1., y_val)
        self.first_drawn.SetMinimum(0.)
        self.first_drawn.SetMaximum(10.)



class ScaleFactorsHisto(varial.tools.FSPlotter):
    pat = p.Suppress('FitChainSumfrom') \
          + p.Word(p.nums).setParseAction(lambda i: float(i[0])) \
          + p.Suppress('to') \
          + p.Word(p.nums).setParseAction(lambda i: float(i[0]))

    def __init__(self, name, legend, with_cuts=False):
        super(ScaleFactorsHisto, self).__init__(name)
        self.with_cuts = with_cuts
        self.legend = legend
        self.canvas_decorators = [IvfEfficiencyLine]

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
            index = leg.index(self.legend)
            res.histo.SetBinContent(i+1, val[index])
            res.histo.SetBinError(i+1, err[index])

        self.result = res
        self.stream_content = [[res]]


summary_chain = varial.tools.ToolChain(
    "Summary", [
        TemplateFitPlots(),
        TemplateFitPlots('TemplateFitPlotsCuts', True),
        TemplatePlots(),
        TemplatePlots('TemplatePlotsCuts', True),
        FitResultCollector(),
        ScaleFactorsHisto('ScaleFactorsHisto', 'B + B Vertex'),
        ScaleFactorsHisto('ScaleFactorsHistoCuts', 'B + B Vertex', True),
        ScaleFactorsHisto('ScaleFactorsHistoBD', 'B + D Vertex'),
        ScaleFactorsHisto('ScaleFactorsHistoBDCuts', 'B + D Vertex', True),
        ScaleFactorsHisto('ScaleFactorsHistoBF', 'B + Fake Vertex'),
        ScaleFactorsHisto('ScaleFactorsHistoBFCuts', 'B + Fake Vertex', True),
    ]
)