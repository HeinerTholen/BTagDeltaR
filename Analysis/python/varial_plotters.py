import itertools
import ROOT
import varial.tools
import varial.generators as gen


def _fix_legend(wrps):
    for w in wrps:
        if 'TTbar' in w.legend:
            w.legend = 'TTbar'
        yield w
jet_plots = varial.tools.FSPlotter(
    'JetPlots',
    filter_keyfunc=lambda w: w.analyzer == 'JetWorker',
    save_lin_log_scale=True,
    hook_loaded_histos=_fix_legend
)


class SampleNormalizer(varial.tools.Tool):
    """Normalize MC cross sections by center of DR distribution. """
    can_reuse = False

    def get_histos_n_factor(self):
        mcee, data = next(gen.fs_mc_stack_n_data_sum(
            lambda w: w.name == 'SelectedPatJetPt' and w.analyzer == 'JetWorker'
        ))
        dh, mh = data.histo, mcee.histo
        #bins = dh.FindBin(-1.51), dh.FindBin(1.51)
        bins = dh.FindBin(20.1), dh.FindBin(100.)
        factor = dh.Integral(*bins) / mh.Integral(*bins)
        canv = next(gen.canvas(
            ((mcee, data),),
            varial.tools.FSPlotter.defaults_attrs['canvas_decorators']
        ))
        return factor, canv

    def run(self):
        # before
        factor, canv = self.get_histos_n_factor()
        next(gen.save_canvas_lin_log([canv], lambda _: 'before'))

        # alter samples
        for s in varial.analysis.mc_samples().itervalues():
            s.lumi /= factor
            s.x_sec /= factor
        for a in varial.analysis.fs_aliases:
            a.lumi /= factor

        # after
        _, canv = self.get_histos_n_factor()
        next(gen.save_canvas_lin_log([canv], lambda _: 'after'))

        self.result = varial.wrappers.FloatWrapper(
            factor,
            name='Lumi factor'
        )


class VtxBeeDeePlotter(varial.tools.FSPlotter):
    def set_up_content(self):
        self.canvas_decorators = [varial.rendering.Legend]
        wrps = self.lookup('../FSHistoLoader')
        wrps = itertools.ifilter(
            lambda w: ('ee' in w.name or 'MatchSig' == w.name)
                      and 'TTbarBDMatch' == w.sample
                      and 'LtDee' not in w.name
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
        'VtxNSharedTracks',
        'LostBeeMatchNDee',
    ] 
    or 'VtxPtLead' in w.name
    or 'VtxPtSubLead' in w.name
)
dist_plotter = varial.tools.FSPlotter(
    'DistPlotter',
    filter_keyfunc=lambda w: (w.name in [
        'VertexMassVsDr',
        'VertexBeeDistLtDeeDist',
        'VertexBeeMassLtDeeMass',
        'VertexBeeDistLtDeeDistOneSigma',
    ] or 'VertexBeeVsDee' in w.name)
    and w.sample == 'TTbarBDMatch'
)
dist_plotter2 = varial.tools.FSPlotter(
    'MassVsDrPlotter',
    filter_keyfunc=lambda w: (w.name in [
        'VertexMassVsDr',
    ]),
)


def _mkchn(analyzer):
    return varial.tools.ToolChainIndie(
        analyzer, [
            varial.tools.FSHistoLoader(None, lambda w: analyzer == w.analyzer),
            stack_plotter,
            beedee_plotter,
            dist_plotter,
            dist_plotter2,
        ]
    )

chains = [
    #_mkchn('IvfMerged'),
    _mkchn('IvfMergedFilt'),
    _mkchn('IvfMergedFiltLt0p2'),
    _mkchn('IvfMergedFiltGt1p0'),
    _mkchn('IvfMergedFiltCuts'),
    _mkchn('IvfB2cMerged'),
    _mkchn('IvfB2cMergedCuts'),
    #_mkchn('IvfB2cMergedFiltCov')
]








