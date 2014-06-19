import itertools
import varial.tools


class VtxBeeDeePlotter(varial.tools.FSPlotter):
    def set_up_content(self):
        wrps = varial.analysis.lookup('../FSHistoLoader')
        wrps = itertools.ifilter(
            lambda w: self.search_string in w.name and 'TTbarBDMatch'==w.sample,
            wrps
        )
        self.stream_content = [wrps]


mass_plotter = VtxBeeDeePlotter(
    'VtxMass',
    input_result_path='../FSHistoLoader',
    search_string='eeMass'
)
num_track_plotter = VtxBeeDeePlotter(
    'VtxNumTracks',
    input_result_path='../FSHistoLoader',
    search_string='eeNumTracks'
)
stack_plotter = varial.tools.FSPlotter(
    input_result_path='../FSHistoLoader'
)

chain_ivf_merged = varial.tools.ToolChain(
    'IvfMerged', [
        varial.tools.FSHistoLoader(None, lambda w: 'IvfMerged' == w.analyzer),
        stack_plotter,
        num_track_plotter,
        mass_plotter
    ]
)
chain_ivf_merged_filt = varial.tools.ToolChain(
    'IvfMergedFilt', [
        varial.tools.FSHistoLoader(None, lambda w: 'IvfMergedFilt'==w.analyzer),
        stack_plotter,
        num_track_plotter,
        mass_plotter
    ]
)
chain_ivf_b2c_merged = varial.tools.ToolChain(
    'IvfB2cMerged', [
        varial.tools.FSHistoLoader(None, lambda w: 'IvfB2cMerged'==w.analyzer),
        stack_plotter,
        num_track_plotter,
        mass_plotter
    ]
)






