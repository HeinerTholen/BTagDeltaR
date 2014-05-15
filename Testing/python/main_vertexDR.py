import ROOT
import itertools
import glob

from varial import fwliteworker, diskio
from MyUtility.PythonUtil.genParticles import final_b_mesons
from DataFormats.FWLite import Events, Handle


# parameters / handles
DR_for_matching = 0.1
h_genParticles = Handle("vector<reco::GenParticle>")
h_pv = Handle("vector<reco::Vertex>")
h_ivf = Handle("vector<reco::Vertex>")


# functions
def my_deltaR(a, b):
    return (
        (a.eta()-b.eta())**2
        + (float(a.phi())-float(b.phi()))**2
    )**.5
DeltaR = ROOT.Math.VectorUtil.DeltaR
deltaR = lambda a, b: DeltaR(a.p4(), b.p4())
deltaR_vec_to_cand = lambda a, b: my_deltaR(a, b.p4())


def get_all_flight_dirs(vertices, primary_vertex):
    def mkfd(vtx):
        pv, sv = primary_vertex, vtx.position()
        return ROOT.Math.XYZVector(
            sv.x() - pv.x(),
            sv.y() - pv.y(),
            sv.z() - pv.z(),
        )
    return list(mkfd(sv) for sv in vertices)


def matching(flightdirs, gen_final_bs, keyfunc, cutvalue):
    combos = itertools.product(flightdirs, gen_final_bs)                 # all combos
    combos = ((fd, bee, keyfunc(fd, bee)) for fd, bee in combos)         # add dR
    combos = list(c for c in combos if c[2] < cutvalue)                  # select < 0.1
    matched_combos = []
    while combos:
        fd, _, _ = combos[0]                                             # start with a flight direction
        bees = list(                                                     # all bees that pair
            b for fd2, b, _ in combos
            if fd == fd2
        )
        cluster = filter(lambda c: c[1] in bees, combos)                 # all combos with these bees
        match = min(cluster, key=lambda c: c[2])                         # select combo with lowest dR
        combos = filter(                                                 # remove all with same fd and b
            lambda c: not (c[0]==match[0] or c[1]==match[1]),
            combos
        )
        matched_combos.append(match)

    return matched_combos


class Worker(fwliteworker.FwliteWorker):
    def __init__(self, name, collection):
        super(Worker, self).__init__(name)
        self.collection = collection

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

        # final B hadron generator particles
        event.getByLabel("genParticles", h_genParticles)
        genParticles = h_genParticles.product()
        fin_bs = final_b_mesons(genParticles)

        # ivf vertices
        event.getByLabel("goodOfflinePrimaryVertices", h_pv)
        event.getByLabel(self.collection, h_ivf)
        flightdirs = get_all_flight_dirs(
            h_ivf.product(), h_pv.product()[0]
        )
        matched_fds = matching(
            flightdirs, fin_bs, deltaR_vec_to_cand, DR_for_matching
        )

        # fill histos for all events
        fs.NumFinalBs.Fill(len(fin_bs))
        fs.NumIvfVertices.Fill(len(flightdirs))

        # fill histos for n_matched == 2
        if len(flightdirs) == 2:
            flightdir_dR = my_deltaR(*flightdirs)
            fs.VertexDR.Fill(flightdir_dR)
            self.vtx_dr_histos[len(matched_fds)].Fill(flightdir_dR)


if __name__ == '__main__':
    workers = [
        Worker("ivf_merged", "inclusiveMergedVertices"),
        Worker("ivf_merged_filt", "inclusiveMergedVerticesFiltered"),
        Worker("ivf_b2c_merged", "bToCharmDecayVertexMerged"),
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
