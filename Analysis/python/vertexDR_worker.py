import glob
import ROOT
ROOT.gROOT.SetBatch()

from varial import fwliteworker, diskio, wrappers
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
    def __init__(self, name, collection, filter_ivf_dr=0.):
        super(Worker, self).__init__(name)
        self.collection = collection
        self.filter_ivf_dr = filter_ivf_dr
        self.result = wrappers.Wrapper(name=name)

    def node_setup(self, init_wrp):
        name = self.name
        if not hasattr(init_wrp, 'announced'):
            init_wrp.announced = True
            print "Starting:", \
                init_wrp.__dict__.get('sample'), init_wrp.filenames

        fs = self.result
        fs.NumIvfVertices = ROOT.TH1D(
            "NumIvfVertices" + "_" + name,
            ";number of IVF vertices;number of events",
            8, -.5, 7.5
        )

        if "Run" in init_wrp.__dict__.get('sample'):
            return

        fs.NumFinalBs = ROOT.TH1D(
            "NumFinalBs" + "_" + name,
            ";number of final B's;number of events",
            8, -.5, 7.5
        )

        # with n_ivf == 2, no / one / two matched
        fs.VertexDR = ROOT.TH1D(
            "VertexDR" + "_" + name,
            ";#Delta R;number of events",
            100, 0., 5.
        )
        fs.DrMomentumFlightdirNoMatch = ROOT.TH1D(
            "DrMomentumFlightdirNoMatch" + "_" + name,
            ";#Delta R;number of vertices",
            100, 0., 1.
        )
        fs.DrMomentumFlightdirOneMatch = ROOT.TH1D(
            "DrMomentumFlightdirOneMatch" + "_" + name,
            ";#Delta R;number of vertices",
            100, 0., 1.
        )
        fs.DrMomentumFlightdirTwoMatch = ROOT.TH1D(
            "DrMomentumFlightdirTwoMatch" + "_" + name,
            ";#Delta R;number of vertices",
            100, 0., 1.
        )
        fs.VertexDRNoMatch = ROOT.TH1D(
            "VertexDRNoMatch" + "_" + name,
            ";#Delta R;number of events",
            100, 0., 5.
        )
        fs.VertexDROneMatch = ROOT.TH1D(
            "VertexDROneMatch" + "_" + name,
            ";#Delta R;number of events",
            100, 0., 5.
        )
        fs.VertexDRTwoMatch = ROOT.TH1D(
            "VertexDRTwoMatch" + "_" + name,
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
        event.getByLabel("goodOfflinePrimaryVertices", h_pv)
        ivf_vtx_fd = get_tuples_with_flight_dirs(
            ivf_vtx, h_pv.product()[0]
        )

        # filter ivf vertices for momentum direction
        if self.filter_ivf_dr:
            ivf_vtx_fd = filter(
                lambda v: deltaR_vec_to_vec(
                    v[1],
                    v[0].p4()
                ) < self.filter_ivf_dr,
                ivf_vtx_fd
            )

        fs.NumIvfVertices.Fill(len(ivf_vtx_fd))

        # end of story for data
        if event.eventAuxiliary().isRealData():
            return

        # final B hadron generator particles
        event.getByLabel("genParticles", h_genParticles)
        genParticles = h_genParticles.product()
        fin_bs = final_b_mesons(genParticles)
        fs.NumFinalBs.Fill(len(fin_bs))

        # right now, use only with exactly two ivf sv's
        if len(ivf_vtx_fd) != 2:
            return

        # if more than two ivf vtx, take the ones closest together
        if len(ivf_vtx_fd) > 2:
            ivf_vtx_fd = min(
                itertools.combinations(ivf_vtx_fd, 2),
                key=lambda c: fd_dr(*c)
            )

        matched = matching(
            ivf_vtx_fd,
            fin_bs,
            lambda a, b: deltaR_vec_to_cand(a[1], b),
            DR_for_matching
        )

        # two or more vertices
        if len(ivf_vtx_fd) > 1:
            n_matched = len(matched)
            dr = fd_dr(*ivf_vtx_fd)
            fs.VertexDR.Fill(dr)
            self.vtx_dr_histos[n_matched].Fill(dr)

            v = ivf_vtx_fd[0]
            self.vtx_dr_mom_fd_histos[n_matched].Fill(
                deltaR_vec_to_vec(v[1], v[0].p4())
            )
            v = ivf_vtx_fd[1]
            self.vtx_dr_mom_fd_histos[n_matched].Fill(
                deltaR_vec_to_vec(v[1], v[0].p4())
            )

    def node_finalize(self, init_wrp):
        if not hasattr(init_wrp, 'announced2'):
            init_wrp.announced2 = True
            print "Done:", \
                init_wrp.__dict__.get('sample'), init_wrp.filenames


workers = [
    Worker("ivf_merged", "inclusiveMergedVertices"),
    Worker("ivf_merged_filt", "inclusiveMergedVerticesFiltered"),
    Worker("ivf_b2c_merged", "bToCharmDecayVertexMerged"),
    #Worker("ivf_merged_momDR", "inclusiveMergedVertices", 0.1),
    #Worker("ivf_merged_filt_momDR", "inclusiveMergedVerticesFiltered", 0.1),
    #Worker("ivf_b2c_merged_momDR", "bToCharmDecayVertexMerged", 0.1),
]


def standalone():
    event_handles = map(
        lambda f: Events("file:%s" % f),
        glob.glob(
            "/nfs/dust/cms/user/tholenhe/samples/"
            "DoubleVtxEffTTdilep/TTdilep_presel_*.root"
        )
    )

    res = fwliteworker.work(workers, event_handles)
    for r in res.values():
        diskio.write(r)


if __name__ == '__main__':
    fwliteworker.work(workers)
