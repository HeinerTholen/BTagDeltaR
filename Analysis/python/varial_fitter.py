import itertools
import ROOT
import varial.tools
import varial.history
import varial.generators as gen


@varial.history.track_history
def get_slice_from_TH2D(wrp, bin_low, bin_high):
    name = wrp.name + 'from%dto%d' % (bin_low, bin_high)
    histo = wrp.histo.ProjectionY(name, bin_low, bin_high)
    histo = histo.Clone()
    histo.SetTitle(wrp.legend)
    return varial.wrappers.HistoWrapper(histo, **wrp.all_info())


def slice_generator(wrps):
    for wrp in wrps:
        yield get_slice_from_TH2D(wrp, 0, 2)
        yield get_slice_from_TH2D(wrp, 10, 80)


class MassHistoSlicer(varial.tools.Tool):
    def run(self):
        wrps = filter(
            lambda w: isinstance(w.histo, ROOT.TH2D),
            self.lookup('../FSHistoLoader')
        )
        self.result = list(slice_generator(wrps))


fitter_chain = varial.tools.ToolChain(
    'VtxMassFitter', [
        varial.tools.FSHistoLoader(
            filter_keyfunc=lambda w: w.name in [
                'VertexBeeMassTemplate',
                'VertexDeeMassTemplate',
                'VertexMassVsDr'
            ]
        ),
        MassHistoSlicer(),
        varial.tools.FSPlotter(
            'TemplatePlots',
            filter_keyfunc=lambda w: 'Template' in w.name
        ),
        varial.tools.FSPlotter(
            'MassSlicePlots',
            input_result_path="../MassHistoSlicer"
        )
    ]
)