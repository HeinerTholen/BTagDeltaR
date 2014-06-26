import itertools
import varial.tools


class VtxBeeDeePlotter(varial.tools.FSPlotter):
    def set_up_content(self):
        self.canvas_decorators = [varial.rendering.Legend]
        wrps = self.lookup('../FSHistoLoader')
        wrps = itertools.ifilter(
            lambda w: ('ee' in w.name or 'MatchSig' == w.name) and 'TTbarBDMatch'==w.sample,
            wrps
        )
        wrps = varial.generators.apply_histo_linecolor(
            wrps, varial.settings.default_colors)
        wrps = sorted(wrps, key=lambda w: w.name[-5:])
        wrps = varial.generators.group(wrps, lambda w: w.name[-5:])
        self.stream_content = wrps


beedee_plotter = VtxBeeDeePlotter()
num_track_plotter = VtxBeeDeePlotter(
    'VtxNumTracks',
    input_result_path='../FSHistoLoader',
    search_string='eeNumTracks'
)
stack_plotter = varial.tools.FSPlotter(
    input_result_path='../FSHistoLoader',
    filter_keyfunc=lambda w: w.name in [
        'VertexDR',
        'NumIvfVertices',
        'DrMomentumFlightdir',
        'VertexDRTwoMatch',
        'NumFinalBs',
    ]
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








