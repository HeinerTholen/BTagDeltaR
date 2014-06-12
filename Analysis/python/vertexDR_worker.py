import glob
import ROOT
ROOT.gROOT.SetBatch()
ROOT.TH1.AddDirectory(False)

from varial import fwliteworker, wrappers
from MyUtility.PythonUtil.genParticles import final_b_mesons
from MyUtility.PythonUtil.eventlooputility import *
from DataFormats.FWLite import Events, Handle


# parameters / handles
DR_for_matching = 0.1
h_genParticles = Handle("vector<reco::GenParticle>")
h_pv = Handle("vector<reco::Vertex>")
h_ivf = Handle("vector<reco::Vertex>")
fd_dr = lambda a, b: deltaR_vec_to_vec(a[1], b[1])


def get_tuples_with_flight_dirs(vertices, primary_vertex):
    """returns [(sv1, fd1), (sv2, fd2), ...]"""
    def mkfd(vtx):
        pv, sv = primary_vertex, vtx.position()
        return ROOT.Math.XYZVector(
            sv.x() - pv.x(),
            sv.y() - pv.y(),
            sv.z() - pv.z(),
        )
    return list((sv, mkfd(sv)) for sv in vertices)


class Worker(fwliteworker.FwliteWorker):
    def __init__(self, name, collection):
        super(Worker, self).__init__(name)
        self.collection = collection
        self.result = wrappers.FileServiceWrapper(name=name)
        self.n_matches_required = -1

    def node_setup(self, init_wrp):
        if not hasattr(init_wrp, 'announced'):
            init_wrp.announced = True
            print "Starting:", init_wrp.sample, \
                init_wrp.event_handle.size(), init_wrp.filenames

        if 'TTbarNoMatch' == init_wrp.sample:
            self.n_matches_required = 0
        elif 'TTbarOneMatch' == init_wrp.sample:
            self.n_matches_required = 1
        elif 'TTbarTwoMatch' == init_wrp.sample:
            self.n_matches_required = 2
        fs = self.result
        fs.NumInEvts = ROOT.TH1D('NumInEvts', 'NumInEvts', 1, 0.5, 1.5)
        fs.NumInEvts.Fill(1., init_wrp.event_handle.size())
        fs.NumIvfVertices = ROOT.TH1D(
            "NumIvfVertices",
            ";number of IVF vertices;number of events",
            8, -.5, 7.5
        )
        fs.VertexDR = ROOT.TH1D(
            "VertexDR",
            ";#Delta R;number of events",
            100, 0., 5.
        )
        fs.DrMomentumFlightdir = ROOT.TH1D(
            "DrMomentumFlightdir",
            ";#Delta R;number of vertices",
            100, 0., 1.
        )

        # from here on only MC histograms
        if "Run" in init_wrp.sample:
            return

        fs.NumFinalBs = ROOT.TH1D(
            "NumFinalBs",
            ";number of final B's;number of events",
            8, -.5, 7.5
        )

        # with n_ivf == 2, no / one / two matched
        fs.DrMomentumFlightdirNoMatch = ROOT.TH1D(
            "DrMomentumFlightdirNoMatch",
            ";#Delta R;number of vertices",
            100, 0., 1.
        )
        fs.DrMomentumFlightdirOneMatch = ROOT.TH1D(
            "DrMomentumFlightdirOneMatch",
            ";#Delta R;number of vertices",
            100, 0., 1.
        )
        fs.DrMomentumFlightdirTwoMatch = ROOT.TH1D(
            "DrMomentumFlightdirTwoMatch",
            ";#Delta R;number of vertices",
            100, 0., 1.
        )
        fs.VertexDRNoMatch = ROOT.TH1D(
            "VertexDRNoMatch",
            ";#Delta R;number of events",
            100, 0., 5.
        )
        fs.VertexDROneMatch = ROOT.TH1D(
            "VertexDROneMatch",
            ";#Delta R;number of events",
            100, 0., 5.
        )
        fs.VertexDRTwoMatch = ROOT.TH1D(
            "VertexDRTwoMatch",
            ";#Delta R;number of events",
            100, 0., 5.
        )
        self.vtx_dr_histos = [
            fs.VertexDRNoMatch,
            fs.VertexDROneMatch,
            fs.VertexDRTwoMatch,
        ]
        self.vtx_dr_mom_fd_histos = [
            fs.DrMomentumFlightdirNoMatch,
            fs.DrMomentumFlightdirOneMatch,
            fs.DrMomentumFlightdirTwoMatch,
        ]

    def node_process_event(self, event):
        fs = self.result

        # ivf vertices
        event.getByLabel(self.collection, h_ivf)
        ivf_vtx = h_ivf.product()

        # flight directions
        event.getByLabel("offlinePrimaryVertices", h_pv)
        ivf_vtx_fd = get_tuples_with_flight_dirs(
            ivf_vtx, h_pv.product()[0]
        )

        # if more than two ivf vtx, take the ones closest together
        ivf_vtx_fd_max2 = ivf_vtx_fd
        if len(ivf_vtx_fd) > 2:
            ivf_vtx_fd_max2 = min(
                itertools.combinations(ivf_vtx_fd, 2),
                key=lambda c: fd_dr(*c)
            )

        # matching
        is_real_data = event.eventAuxiliary().isRealData()
        matched = []
        fin_bs = []
        if not is_real_data:
            event.getByLabel("genParticles", h_genParticles)
            genParticles = h_genParticles.product()
            fin_bs = final_b_mesons(genParticles)
            matched = matching(
                ivf_vtx_fd_max2,
                fin_bs,
                lambda a, b: deltaR_vec_to_cand(a[1], b),
                DR_for_matching
            )
        n_matched = len(matched)

        # check for right ttbar or different sample
        if self.n_matches_required not in (-1, n_matched):
            return

        # fill histograms
        fs.NumIvfVertices.Fill(len(ivf_vtx_fd))
        for vtx, fd in ivf_vtx_fd:
            fs.DrMomentumFlightdir.Fill(deltaR_vec_to_vec(fd, vtx.p4()))
        if not is_real_data:
            fs.NumFinalBs.Fill(len(fin_bs))
            for a, b in ivf_vtx_fd:
                self.vtx_dr_mom_fd_histos[n_matched].Fill(
                    deltaR_vec_to_vec(a.p4(), b)
                )

        # dr
        if len(ivf_vtx_fd_max2) == 2:
            dr = fd_dr(*ivf_vtx_fd_max2)
            fs.VertexDR.Fill(dr)
            if not is_real_data:
                self.vtx_dr_histos[n_matched].Fill(dr)

    def node_finalize(self, init_wrp):
        if not hasattr(init_wrp, 'announced2'):
            init_wrp.announced2 = True
            print "Done:", init_wrp.sample, \
                init_wrp.event_handle.size(), init_wrp.filenames


workers = [
    Worker("IvfMerged", "inclusiveMergedVertices"),
    Worker("IvfMergedFilt", "inclusiveMergedVerticesFiltered"),
    Worker("IvfB2cMerged", "bToCharmDecayVertexMerged"),
]


if __name__ == '__main__':
    fwliteworker.work(workers)
