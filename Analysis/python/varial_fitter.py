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
        for i in xrange(0, 100, 20):
            yield get_slice_from_TH2D(wrp, i, i+20)


class MassHistoSlicer(varial.tools.Tool):
    def run(self):
        wrps = list(self.lookup('../FSHistoLoader'))
        self.result = list(slice_generator(wrps))


fitter_chain = varial.tools.ToolChain(
    'VtxMassFitter', [
        varial.tools.FSHistoLoader(
            None,
            lambda w: 'VertexMassVsDr' == w.name
        ),
        MassHistoSlicer(),
        varial.tools.FSPlotter(
            'MassSlicePlots',
            input_result_path="../MassHistoSlicer"
        )
    ]
)