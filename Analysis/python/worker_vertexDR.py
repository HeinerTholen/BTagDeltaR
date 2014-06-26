import ROOT
ROOT.gROOT.SetBatch()
ROOT.TH1.AddDirectory(False)

from varial import fwliteworker, wrappers
from MyUtility.PythonUtil.genParticles import *
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


class PreWorker(fwliteworker.FwliteWorker):
    def node_process_event(self, event):
        is_real_data = event.eventAuxiliary().isRealData()
        if not is_real_data:
            event.getByLabel("genParticles", h_genParticles)
            event.gen_particles = h_genParticles.product()
            event.fin_b = final_b_hadrons(event.gen_particles)


class Worker(fwliteworker.FwliteWorker):
    def __init__(self, name, collection, filter_vtx=False, dee_cov_match=False):
        super(Worker, self).__init__(name)
        self.collection = collection
        self.n_matches_required = -1
        self.filter_vtx = filter_vtx
        self.dee_cov_match = dee_cov_match

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
        elif 'TTbarBDMatch' == init_wrp.sample:
            self.n_matches_required = 3
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
        fs.make(
            "VtxPtLeadMass",
            ";pt leading vertex: mass; number of events",
            100, 0., 10.
        )
        fs.make(
            "VtxPtSubLeadMass",
            ";pt sub leading vertex: mass; number of events",
            100, 0., 10.
        )
        fs.make(
            "VtxPtLeadPt",
            ";pt leading vertex: p_T; number of events",
            100, 0., 500.
        )
        fs.make(
            "VtxPtSubLeadPt",
            ";pt sub leading vertex: p_T; number of events",
            100, 0., 500.
        )
        fs.make(
            "VtxPtLeadEta",
            ";pt leading vertex: #eta; number of events",
            100, -5., 5.
        )
        fs.make(
            "VtxPtSubLeadEta",
            ";pt sub leading vertex: #eta; number of events",
            100, -5., 5.
        )
        fs.make(
            "VtxPtLeadNumTracks",
            ";pt leading vertex: no. tracks; number of events",
            11, -.5, 10.5
        )
        fs.make(
            "VtxPtSubLeadNumTracks",
            ";pt sub leading vertex: no. tracks; number of events",
            11, -.5, 10.5
        )
 
        # from here on only MC histograms
        if "Run" in init_wrp.sample:
            return

        fs.NumFinalBs = ROOT.TH1D(
            "NumFinalBs",
            ";number of final B's;number of events",
            8, -.5, 7.5
        )
        fs.make(
            "DrFdFinalBs",
            ";#Delta R; number of events",
            100, 0., 1.
        )
        fs.make(
            "DrMomFinalBs",
            ";#Delta R; number of events",
            100, 0., 1.
        )

        # for B / D vertices
        fs.make(
            "VtxBeeNumTracks",
            ";number of Tracks;number of vertices",
            16, -.5, 15.5
        )
        fs.make(
            "VtxDeeNumTracks",
            ";number of Tracks;number of vertices",
            16, -.5, 15.5
        )
        fs.make(
            "VtxBeeMass",
            ";vertex mass;number of vertices",
            100, 0., 10.
        )
        fs.make(
            "VtxDeeMass",
            ";vertex mass;number of vertices",
            100, 0., 10.
        )
        fs.make(
            "VtxBeeDeeMatchSig",
            ";#DeltaX for matching B and D;number of vertices",
            100, 0., 100
        )

        # DrMomentumFlightdir
        #fs.DrMomentumFlightdirNoMatch = ROOT.TH1D(
        #    "DrMomentumFlightdirNoMatch",
        #    ";#Delta R;number of vertices",
        #    100, 0., 1.
        #)
        #fs.DrMomentumFlightdirOneMatch = ROOT.TH1D(
        #    "DrMomentumFlightdirOneMatch",
        #    ";#Delta R;number of vertices",
        #    100, 0., 1.
        #)
        #fs.DrMomentumFlightdirBDMatch = ROOT.TH1D(
        #    "DrMomentumFlightdirOneMatch",
        #    ";#Delta R;number of vertices",
        #    100, 0., 1.
        #)
        #fs.DrMomentumFlightdirTwoMatch = ROOT.TH1D(
        #    "DrMomentumFlightdirTwoMatch",
        #    ";#Delta R;number of vertices",
        #    100, 0., 1.
        #)

        ## VertexDR
        #fs.VertexDRNoMatch = ROOT.TH1D(
        #    "VertexDRNoMatch",
        #    ";#Delta R;number of events",
        #    100, 0., 5.
        #)
        #fs.VertexDROneMatch = ROOT.TH1D(
        #    "VertexDROneMatch",
        #    ";#Delta R;number of events",
        #    100, 0., 5.
        #)
        #fs.VertexDRBDMatch = ROOT.TH1D(
        #    "VertexDRBDMatch",
        #    ";#Delta R;number of events",
        #    100, 0., 5.
        #)
        #fs.VertexDRTwoMatch = ROOT.TH1D(
        #    "VertexDRTwoMatch",
        #    ";#Delta R;number of events",
        #    100, 0., 5.
        #)
        #self.vtx_dr_histos = [
        #    fs.VertexDRNoMatch,
        #    fs.VertexDROneMatch,
        #    fs.VertexDRTwoMatch,
        #    fs.VertexDRBDMatch,
        #]
        #self.vtx_dr_mom_fd_histos = [
        #    fs.DrMomentumFlightdirNoMatch,
        #    fs.DrMomentumFlightdirOneMatch,
        #    fs.DrMomentumFlightdirTwoMatch,
        #    fs.DrMomentumFlightdirBDMatch,
        #]

    def node_process_event(self, event):
        fs = self.result

        # ivf vertices
        event.getByLabel(self.collection, h_ivf)
        ivf_vtx = h_ivf.product()
        if self.filter_vtx:
            ivf_vtx = filter(
                lambda v: (
                    abs(v.p4().eta()) < 2
                    and v.p4().pt() > 8
                    and 6.5 > v.p4().mass() > 1.4
                    and v.nTracks() > 2
                ),
                ivf_vtx
            )
        ivf_vtx_pt_sort = sorted(ivf_vtx, key=lambda v: -v.p4().pt())

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
        matched_bd = []
        fin_b = []
        fin_bd = []
        if not is_real_data:
            gen_particles = event.gen_particles
            fin_b = event.fin_b

            matched = matching(
                ivf_vtx_fd_max2,
                fin_b,
                #lambda a, b: deltaR_vec_to_cand(a[1], b),
                lambda a, b: deltaR_cand_to_cand(a[0], b),
                DR_for_matching
            )

            if len(matched) == 1:
                # try to match D hadron as well
                fin_b_match = matched[0][1]
                fin_d = final_d_hadrons(get_all_daughters(
                    gen_particles, [fin_b_match]))
                fin_bd = [fin_b_match] + fin_d
                if self.dee_cov_match:
                    matched_bd = matching(
                        ivf_vtx_fd_max2,
                        fin_bd,
                        lambda a, b: covariance_significance(a[0], b),
                        10.
                    )
                else:
                    matched_bd = matching(
                        ivf_vtx_fd_max2,
                        fin_bd,
                        lambda a, b: deltaR_cand_to_cand(a[0], b),
                        DR_for_matching
                    )

                # discard, if b not present anymore, else sort b to front
                if fin_b_match not in (gp for _, gp, _ in matched_bd):
                    matched_bd = []
                elif 2 == len(matched_bd):
                    if matched_bd[0][1] != fin_b_match:
                        dee, bee = matched_bd
                        matched_bd = bee, dee

        n_matched = len(matched)
        if 1 == n_matched and 2 == len(matched_bd):
            n_matched = 3

        # check for right ttbar or different sample
        if self.n_matches_required not in (-1, n_matched):
            return

        # fill histograms
        fs.NumIvfVertices.Fill(len(ivf_vtx_fd))
        if ivf_vtx_pt_sort:
            fs.VtxPtLeadMass.Fill(ivf_vtx_pt_sort[0].p4().mass())
            fs.VtxPtLeadPt.Fill(ivf_vtx_pt_sort[0].p4().pt())
            fs.VtxPtLeadEta.Fill(ivf_vtx_pt_sort[0].p4().eta())
            fs.VtxPtLeadNumTracks.Fill(ivf_vtx_pt_sort[0].nTracks())
        if len(ivf_vtx_pt_sort) > 1:
            fs.VtxPtSubLeadMass.Fill(ivf_vtx_pt_sort[1].p4().mass())
            fs.VtxPtSubLeadPt.Fill(ivf_vtx_pt_sort[1].p4().pt())
            fs.VtxPtSubLeadEta.Fill(ivf_vtx_pt_sort[1].p4().eta())
            fs.VtxPtSubLeadNumTracks.Fill(ivf_vtx_pt_sort[1].nTracks())

        for vtx, fd in ivf_vtx_fd:
            fs.DrMomentumFlightdir.Fill(deltaR_vec_to_vec(fd, vtx.p4()))
        if not is_real_data:
            fs.NumFinalBs.Fill(len(fin_b))
            if len(fin_b) > 1:
                fs.DrFdFinalBs.Fill(min(
                    deltaR_vec_to_vec(
                        mkrtvec(mkvec(a.daughter(0).vertex()) - mkvec(a.vertex())), 
                        mkrtvec(mkvec(b.daughter(0).vertex()) - mkvec(b.vertex()))
                    )
                    for a, b in itertools.combinations(fin_b, 2)
                ))
                fs.DrMomFinalBs.Fill(min(
                    deltaR_cand_to_cand(a, b)
                    for a, b in itertools.combinations(fin_b, 2)
                ))
            #for a, b in ivf_vtx_fd:
            #    self.vtx_dr_mom_fd_histos[n_matched].Fill(
            #        deltaR_vec_to_vec(a.p4(), b)
            #    )

        # dr (needs two vertices)
        if len(ivf_vtx_fd_max2) == 2:
            dr = fd_dr(*ivf_vtx_fd_max2)
            fs.VertexDR.Fill(dr)
            #if not is_real_data:
            #    self.vtx_dr_histos[n_matched].Fill(dr)

        # fill B / D vertex histos if in bd mode:
        if 3 == n_matched:
            bee, dee = matched_bd
            bee = bee[0][0]
            dee = dee[0][0]
            fs.VtxBeeNumTracks.Fill(bee.nTracks())
            fs.VtxDeeNumTracks.Fill(dee.nTracks())
            fs.VtxBeeMass.Fill(bee.p4().M())
            fs.VtxDeeMass.Fill(dee.p4().M())

            # matching significance
            matching_bd_cov = matching(
                ivf_vtx_fd_max2,
                fin_bd,
                lambda a, b: covariance_significance(a[0], b),
                100.
            )
            for _, _, sig in matching_bd_cov:
                fs.VtxBeeDeeMatchSig.Fill(sig)

    def node_finalize(self, init_wrp):
        if not hasattr(init_wrp, 'announced2'):
            init_wrp.announced2 = True
            print "Done:", init_wrp.sample, \
                init_wrp.event_handle.size(), init_wrp.filenames


workers = [
    PreWorker("PreWorker"),
    #Worker("IvfMerged", "inclusiveMergedVertices"),
    Worker("IvfMergedFilt", "inclusiveMergedVerticesFiltered"),
    Worker("IvfB2cMerged", "bToCharmDecayVertexMerged"),
    Worker("IvfB2cMergedFilt", "bToCharmDecayVertexMergedFilt", True),
    #Worker("IvfB2cMergedFiltCov", "bToCharmDecayVertexMergedFilt", True, True),
]


if __name__ == '__main__':
    fwliteworker.work(workers)
