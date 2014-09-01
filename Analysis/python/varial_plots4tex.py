import os
from os.path import join

import varial.tools


control_plot_names = [
    'VtxPtLeadEta',
    'VtxPtSubLeadEta',
    'VtxPtLeadPt',
    'VtxPtSubLeadPt',

    'VtxPtLeadMass',
    'VtxPtSubLeadMass',
    'VtxPtLeadNumTracks',
    'VtxPtSubLeadNumTracks',

    'VertexFdDR',
    'VertexMomDR',
    'NumIvfVertices',
    'DrMomentumFlightdir',
]
control_plot_collections = [
    'IvfMergedFilt',
    'IvfB2cMerged',
    'IvfB2cMergedCuts',
]
include_str = r"\includegraphics[width=0.49\textwidth]{%s}"


class ControlPlotConverter(varial.tools.Tool):
    def __init__(self, collection):
        super(ControlPlotConverter, self).__init__(collection)

    def run(self):
        path = varial.analysis.lookup_path(
            '../../' + self.name + '/ControlPlots')
        for plot in control_plot_names:
            p = join(path, self.name + '_' + plot)
            os.system('convert %s.eps %s.pdf' % (p, p))
            os.system('cp %s.pdf %s' % (p, self.cwd))


class TexCrtlPlts(varial.tools.ToolChain):
    def __init__(self):
        super(TexCrtlPlts, self).__init__(None, list(
            ControlPlotConverter(c) for c in control_plot_collections
        ))

    def run(self):
        super(TexCrtlPlts, self).run()
        N = 4
        for c in control_plot_collections:
            for i_group, group in enumerate(
                    control_plot_names[k:k+N]
                    for k in range(0, len(control_plot_names), N)
            ):
                include_strs = list(
                    include_str % ("%s/%s/%s_%s" % (self.name, c, c, name)) + '\n'
                    for name in group
                )
                fname = join(varial.analysis.cwd, c) + "_%i.tex" % i_group
                with open(fname, 'w') as f:
                    f.writelines(include_strs)


