import ROOT
ROOT.gROOT.SetBatch()
ROOT.TH1.AddDirectory(False)

from libMyUtilityPythonUtil import get_2d_3d_distances
from varial import fwliteworker, wrappers
from MyUtility.PythonUtil.genParticles import *
from MyUtility.PythonUtil.eventlooputility import *
from DataFormats.FWLite import Events, Handle


# parameters / handles
DR_for_matching = 0.1
h_gen_particles = Handle("vector<reco::GenParticle>")
h_pu_weight = Handle("double")
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
        event.weight = 1.
        if not is_real_data:
            event.getByLabel("genParticles", h_gen_particles)
            event.gen_particles = h_gen_particles.product()
            event.fin_b = final_b_hadrons(event.gen_particles)

            event.getByLabel("puWeight", "PUWeightTrue", h_pu_weight)
            event.weight *= h_pu_weight.product()[0]


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
        fs.make(
            "VertexFdDR",
            ";#Delta R;number of events",
            100, 0., 5.
        )
        fs.make(
            "VertexMomDR",
            ";#Delta R;number of events",
            100, 0., 5.
        )
        fs.make(
            "DrMomentumFlightdir",
            ";#Delta R;number of vertices",
            100, 0., 1.
        )
        fs.VertexMassVsDr = ROOT.TH2D(
            'VertexMassVsDr',
            ';Vertex #Delta R; Vertex Mass',
            100, 0., 5.,
            100, 0., 10.
        )
        fs.make(
            'VertexBeeMassTemplate',
            ";B vertex candidate mass; number of events",
            100, 0., 10.
        )
        fs.make(
            'VertexDeeMassTemplate',
            ";D vertex candidate mass; number of events",
            100, 0., 10.
        )

        # control plots
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
            ";pt leading vertex: p_{T}; number of events",
            100, 0., 500.
        )
        fs.make(
            "VtxPtSubLeadPt",
            ";pt sub leading vertex: p_{T}; number of events",
            100, 0., 500.
        )
        fs.make(
            "VtxPtLeadEta",
            ";pt leading vertex: momentum #eta; number of events",
            100, -2.5, 2.5
        )
        fs.make(
            "VtxPtSubLeadEta",
            ";pt sub leading vertex: momentum #eta; number of events",
            100, -2.5, 2.5
        )
        fs.make(
            "VtxPtLeadFdEta",
            ";pt leading vertex: flight direction #eta; number of events",
            100, -2.5, 2.5
        )
        fs.make(
            "VtxPtSubLeadFdEta",
            ";pt sub leading vertex: flight direction #eta; number of events",
            100, -2.5, 2.5
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

        fs.make(
            "EventWeight",
            ";MC event weight;number of events",
            100, -5., 5.
        )
        fs.make(
            "NumFinalBs",
            ";number of final B's;number of events",
            8, -.5, 7.5
        )
        fs.make(
            "DrFdFinalBs",
            ";#Delta R; number of events",
            100, 0., 5.
        )
        fs.make(
            "DrMomFinalBs",
            ";#Delta R; number of events",
            100, 0., 5.
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
        fs.VertexBeeDistVsDeeDist = ROOT.TH2D(
            'VertexBeeDistVsDeeDist',
            ';D Vertex Dist3D; B Vertex Dist3D',
            100, 0., 5.,
            100, 0., 5.
        )

    def node_process_event(self, event):
        fs = self.result

        # ivf vertices
        event.getByLabel(self.collection, h_ivf)
        ivf_vtx = list(h_ivf.product())

        # pv
        event.getByLabel("offlinePrimaryVertices", h_pv)
        pv = h_pv.product()[0]

        # flight distances
        for sv in ivf_vtx:
            sv.__dict__.update(get_2d_3d_distances(sv, pv))

        # flight directions
        ivf_vtx_fd = get_tuples_with_flight_dirs(
            ivf_vtx, pv
        )
        if self.filter_vtx:
            ivf_vtx_fd = filter(
                lambda v: (
                    abs(v[0].p4().eta()) < 2
                    and v[0].p4().pt() > 8
                    and 6.5 > v[0].p4().mass() > 1.4
                    and v[0].nTracks() > 2
                    and deltaR_vec_to_vec(v[1], v[0].p4()) < 0.1
                    and v[0].dxy_val < 2.5
                    and (v[0].dxy_val / v[0].dxy_err) > 3.
                    and (v[0].dxyz_val / v[0].dxyz_err) > 5.
                ),
                ivf_vtx_fd
            )

        ivf_vtx_fd_pt_sort = sorted(ivf_vtx_fd, key=lambda v: -v[0].p4().pt())

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

        ################################################### fill histograms ###
        w = event.weight
        fs.NumIvfVertices.Fill(len(ivf_vtx_fd), w)

        if ivf_vtx_fd_pt_sort:
            fs.VtxPtLeadMass.Fill(ivf_vtx_fd_pt_sort[0][0].p4().mass(), w)
            fs.VtxPtLeadPt.Fill(ivf_vtx_fd_pt_sort[0][0].p4().pt(), w)
            fs.VtxPtLeadEta.Fill(ivf_vtx_fd_pt_sort[0][0].p4().eta(), w)
            fs.VtxPtLeadFdEta.Fill(ivf_vtx_fd_pt_sort[0][1].eta(), w)
            fs.VtxPtLeadNumTracks.Fill(ivf_vtx_fd_pt_sort[0][0].nTracks(), w)
        if len(ivf_vtx_fd_pt_sort) > 1:
            fs.VtxPtSubLeadMass.Fill(ivf_vtx_fd_pt_sort[1][0].p4().mass(), w)
            fs.VtxPtSubLeadPt.Fill(ivf_vtx_fd_pt_sort[1][0].p4().pt(), w)
            fs.VtxPtSubLeadEta.Fill(ivf_vtx_fd_pt_sort[1][0].p4().eta(), w)
            fs.VtxPtSubLeadFdEta.Fill(ivf_vtx_fd_pt_sort[1][1].eta(), w)
            fs.VtxPtSubLeadNumTracks.Fill(ivf_vtx_fd_pt_sort[1][0].nTracks(), w)

        for vtx, fd in ivf_vtx_fd:
            fs.DrMomentumFlightdir.Fill(deltaR_vec_to_vec(fd, vtx.p4()), w)

        # dr (needs two vertices)
        if len(ivf_vtx_fd_max2) == 2:
            dr_fd = fd_dr(*ivf_vtx_fd_max2)
            dr_mom = deltaR_cand_to_cand(ivf_vtx_fd_max2[0][0],
                                         ivf_vtx_fd_max2[1][0])
            fs.VertexFdDR.Fill(dr_fd, w)
            fs.VertexMomDR.Fill(dr_mom, w)

            # source for fitted histograms
            fs.VertexMassVsDr.Fill(
                dr_mom,
                ivf_vtx_fd_max2[0][0].p4().mass(),
                w
            )
            fs.VertexMassVsDr.Fill(
                dr_mom,
                ivf_vtx_fd_max2[1][0].p4().mass(),
                w
            )

            # fit templates
            if dr_mom < 0.05:
                ivf_dist_sorted = sorted(ivf_vtx_fd_max2,
                                         key=lambda v: v[0].dxyz_val)
                fs.VertexBeeMassTemplate.Fill(
                    ivf_dist_sorted[0][0].p4().mass(), w)
                fs.VertexDeeMassTemplate.Fill(
                    ivf_dist_sorted[1][0].p4().mass(), w)

        ######################################################### MC histos ###
        if not is_real_data:
            fs.EventWeight.Fill(w)
            fin_b = filter(lambda b: b.p4().pt() > 15.
                                     and abs(b.p4().eta()) < 2.1, fin_b)
            fs.NumFinalBs.Fill(len(fin_b), w)
            if len(fin_b) > 1:
                fs.DrFdFinalBs.Fill(min(
                    deltaR_vec_to_vec(
                        mkrtvec(
                            mkvec(a.daughter(0).vertex()) - mkvec(a.vertex())),
                        mkrtvec(
                            mkvec(b.daughter(0).vertex()) - mkvec(b.vertex()))
                    )
                    for a, b in itertools.combinations(fin_b, 2)
                ))
                fs.DrMomFinalBs.Fill(min(
                    deltaR_cand_to_cand(a, b)
                    for a, b in itertools.combinations(fin_b, 2)
                ))

        # fill B / D vertex histos if in bd mode:
        if 3 == n_matched:
            bee, dee = matched_bd
            bee = bee[0][0]
            dee = dee[0][0]
            fs.VtxBeeNumTracks.Fill(bee.nTracks(), w)
            fs.VtxDeeNumTracks.Fill(dee.nTracks(), w)
            fs.VtxBeeMass.Fill(bee.p4().M(), w)
            fs.VtxDeeMass.Fill(dee.p4().M(), w)
            fs.VertexBeeDistVsDeeDist.Fill(dee.dxyz_val, bee.dxyz_val, w)
            # matching significance
            matching_bd_cov = matching(
                ivf_vtx_fd_max2,
                fin_bd,
                lambda a, b: covariance_significance(a[0], b),
                100.
            )
            for _, _, sig in matching_bd_cov:
                fs.VtxBeeDeeMatchSig.Fill(sig, w)

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
    Worker("IvfB2cMergedFilt", "bToCharmDecayVertexMerged", True),
    #Worker("IvfB2cMergedFiltCov", "bToCharmDecayVertexMergedFilt", True, True),
]


if __name__ == '__main__':
    fwliteworker.work(workers)


