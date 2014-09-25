import itertools
import ROOT
import varial.tools
import varial.generators as gen


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
            varial.tools.FSHistoLoader(
                None, lambda w: analyzer == w.analyzer, io=varial.diskio),
            stack_plotter,
            beedee_plotter,
            dist_plotter,
            dist_plotter2,
        ]
    )

chains = [
    #_mkchn('IvfMerged'),
    _mkchn('IvfMergedFilt'),
    #_mkchn('IvfMergedFiltLt0p2'),
    #_mkchn('IvfMergedFiltGt1p0'),
    _mkchn('IvfMergedFiltCuts'),
    _mkchn('IvfB2cMerged'),
    _mkchn('IvfB2cMergedCuts'),
    #_mkchn('IvfB2cMergedFiltCov')
]








