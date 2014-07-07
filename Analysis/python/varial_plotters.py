import itertools
import ROOT
import varial.tools
import varial.generators as gen


class DaNormalizer(varial.tools.Tool):
    """Normalize MC cross sections by center of DR distribution. """
    def get_histos_n_factor(self):
        mcee, data = next(gen.fs_mc_stack_n_data_sum(
            lambda w: w.name=='VertexMomDR' and w.analyzer=='IvfB2cMergedFilt'
        ))
        dh, mh = data.histo, mcee.histo
        bins = dh.FindBin(1.51), dh.FindBin(3.51)
        factor = dh.Integral(*bins) / mh.Integral(*bins)
        canv = next(gen.canvas(
            ((mcee, data),),
            varial.tools.FSPlotter.defaults_attrs['canvas_decorators']
        ))
        return factor, canv

    def run(self):

        # before
        factor, canv = self.get_histos_n_factor()
        next(gen.save([canv], lambda _: 'before'))

        # alter samples
        for s in varial.analysis.mc_samples().itervalues():
            s.lumi *= factor
            s.x_sec *= factor

        # after
        _, canv = self.get_histos_n_factor()
        next(gen.save([canv], lambda _: 'after'))


class VtxBeeDeePlotter(varial.tools.FSPlotter):
    def set_up_content(self):
        self.canvas_decorators = [varial.rendering.Legend]
        wrps = self.lookup('../FSHistoLoader')
        wrps = itertools.ifilter(
            lambda w: ('ee' in w.name or 'MatchSig' == w.name)
                      and 'TTbarBDMatch' == w.sample
                      and type(w.histo) == ROOT.TH1D,
            wrps
        )
        wrps = gen.apply_histo_linecolor(
            wrps, varial.settings.default_colors)
        wrps = sorted(wrps, key=lambda w: w.name[-5:])
        wrps = gen.group(wrps, lambda w: w.name[-5:])
        self.stream_content = wrps
beedee_plotter = VtxBeeDeePlotter()


stack_plotter = varial.tools.FSPlotter(
    "ControlPlots",
    input_result_path='../FSHistoLoader',
    filter_keyfunc=lambda w: w.name in [
        'VertexMomDR',
        'VertexFdDR',
        'NumIvfVertices',
        'DrMomentumFlightdir',
        'VertexDRTwoMatch',
        'NumFinalBs',
        'DrFdFinalBs',
        'DrMomFinalBs',
        'VertexDrFdMomDiff',
        'EventWeight',
    ] 
    or 'VtxPtLead' in w.name
    or 'VtxPtSubLead' in w.name
)

chain_ivf_merged = varial.tools.ToolChain(
    'IvfMerged', [
        varial.tools.FSHistoLoader(None, lambda w: 'IvfMerged' == w.analyzer),
        stack_plotter,
        beedee_plotter,
    ]
)
chain_ivf_merged_filt = varial.tools.ToolChain(
    'IvfMergedFilt', [
        varial.tools.FSHistoLoader(None, lambda w: 'IvfMergedFilt'==w.analyzer),
        stack_plotter,
        beedee_plotter,
    ]
)
chain_ivf_b2c_merged = varial.tools.ToolChain(
    'IvfB2cMerged', [
        varial.tools.FSHistoLoader(None, lambda w: 'IvfB2cMerged'==w.analyzer),
        stack_plotter,
        beedee_plotter,
    ]
)
chain_ivf_b2c_merged_filt = varial.tools.ToolChain(
    'IvfB2cMergedFilt', [
        varial.tools.FSHistoLoader(None, lambda w: 'IvfB2cMergedFilt'==w.analyzer),
        stack_plotter,
        beedee_plotter,
    ]
)
chain_ivf_b2c_merged_filt_cov = varial.tools.ToolChain(
    'IvfB2cMergedFiltCov', [
        varial.tools.FSHistoLoader(None, lambda w: 'IvfB2cMergedFiltCov'==w.analyzer),
        stack_plotter,
        beedee_plotter,
    ]
)








