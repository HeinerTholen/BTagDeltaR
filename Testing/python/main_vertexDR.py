import ROOT
import glob

from varial import fwliteworker, diskio
from MyUtility.PythonUtil.genParticles import final_b_mesons
from MyUtility.PythonUtil.eventlooputility import *
from DataFormats.FWLite import Events, Handle


# parameters / handles
DR_for_matching = 0.1
h_genParticles = Handle("vector<reco::GenParticle>")
h_pv = Handle("vector<reco::Vertex>")
h_ivf = Handle("vector<reco::Vertex>")


def get_all_flight_dirs(vertices, primary_vertex):
    def mkfd(vtx):
        pv, sv = primary_vertex, vtx.position()
        return ROOT.Math.XYZVector(
            sv.x() - pv.x(),
            sv.y() - pv.y(),
            sv.z() - pv.z(),
        )
    return list(mkfd(sv) for sv in vertices)


class Worker(fwliteworker.FwliteWorker):
    def __init__(self, name, collection, filter_ivf_dr=0.):
        super(Worker, self).__init__(name)
        self.collection = collection
        self.filter_ivf_dr = filter_ivf_dr

        fs = self.result
        fs.NumFinalBs = ROOT.TH1D(
            "NumFinalBs" + "_" + name,
            ";number of final B's;number of events",
            8, -.5, 7.5
        )
        fs.NumIvfVertices = ROOT.TH1D(
            "NumIvfVertices" + "_" + name,
            ";number of IVF vertices;number of events",
            8, -.5, 7.5
        )

        # with n_ivf == 2, no / one / two matched
        fs.VertexDR = ROOT.TH1D(
            "VertexDR" + "_" + name,
            ";#Delta R;number of events",
            100, 0., 5.
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
            fs.VertexDRTwoMatch
        ]

    def node_setup(self, event_handle):
        print "Starting:", event_handle._filenames, self.result.name

    def node_process_event(self, event):
        fs = self.result

        # ivf vertices
        event.getByLabel(self.collection, h_ivf)
        ivf_vtx = h_ivf.product()
        if self.filter_ivf_dr:
            ivf_vtx = filter(
                lambda v: my_deltaR(v.position(), v.p4()) < self.filter_ivf_dr,
                ivf_vtx
            )

        # final B hadron generator particles
        event.getByLabel("genParticles", h_genParticles)
        genParticles = h_genParticles.product()
        fin_bs = final_b_mesons(genParticles)

        # flight directions
        event.getByLabel("goodOfflinePrimaryVertices", h_pv)

        flightdirs = get_all_flight_dirs(
            ivf_vtx, h_pv.product()[0]
        )
        matched_fds = matching(
            flightdirs, fin_bs, deltaR_vec_to_cand, DR_for_matching
        )

        # fill histos for all events
        fs.NumFinalBs.Fill(len(fin_bs))
        fs.NumIvfVertices.Fill(len(flightdirs))

        # two or more vertices
        if len(flightdirs) > 1:

            # if more than two, take the ones closer together
            if len(flightdirs) > 2:
                combos = itertools.combinations(flightdirs, 2)
                flightdirs = min(combos, lambda a, b: my_deltaR(a, b))

            flightdir_dr = my_deltaR(*flightdirs)
            fs.VertexDR.Fill(flightdir_dr)
            self.vtx_dr_histos[len(matched_fds)].Fill(flightdir_dr)


if __name__ == '__main__':
    workers = [
        Worker("ivf_merged", "inclusiveMergedVertices"),
        Worker("ivf_merged_filt", "inclusiveMergedVerticesFiltered"),
        Worker("ivf_b2c_merged", "bToCharmDecayVertexMerged"),
        Worker("ivf_merged", "inclusiveMergedVertices", 0.1),
        Worker("ivf_merged_filt", "inclusiveMergedVerticesFiltered", 0.1),
        Worker("ivf_b2c_merged", "bToCharmDecayVertexMerged", 0.1),
    ]

    event_handles = map(
        lambda f: Events("file:%s" % f),
        glob.glob(
            "/nfs/dust/cms/user/tholenhe/samples/"
            "DoubleVtxEffTTdilep/TTdilep_presel_*.root"
        )[:10]
    )

    res = fwliteworker.work(workers, event_handles)
    for r in res.values():
        diskio.write(r)
