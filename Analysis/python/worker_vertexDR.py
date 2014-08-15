import ROOT
ROOT.gROOT.SetBatch()
ROOT.TH1.AddDirectory(False)

from libMyUtilityPythonUtil import get_2d_3d_distances
from varial import fwliteworker
from MyUtility.PythonUtil.genParticles import *
from MyUtility.PythonUtil.eventlooputility import *
from DataFormats.FWLite import Events, Handle


# globals
DR_for_matching = 0.1
fd_dr = lambda a, b: deltaR_vec_to_vec(a[1], b[1])
h_gen_particles = Handle("vector<reco::GenParticle>")
h_pu_weight = Handle("double")
h_pv = Handle("vector<reco::Vertex>")
h_ivf_normal = Handle("vector<reco::Vertex>")
h_ivf_merged = Handle("vector<reco::MergedVertex>")
h_jets = Handle("vector<pat::Jet>")


def get_tuples_with_flight_dirs(vertices, primary_vertex):
    """returns [(sv1, fd1), (sv2, fd2), ...]"""
    def mkfd(vtx):
        pv, sv = primary_vertex.position(), vtx.position()
        return ROOT.Math.XYZVector(
            sv.x() - pv.x(),
            sv.y() - pv.y(),
            sv.z() - pv.z(),
        )
    return list((sv, mkfd(sv)) for sv in vertices)


def deltaR_from_pv(sv_fd, gen_part, gen_part_at_pv):
    sv = gen_part.daughter(0)
    pv = gen_part_at_pv.mother()
    gen_fd = ROOT.Math.XYZVector(
        sv.vx() - pv.vx(),
        sv.vy() - pv.vy(),
        sv.vz() - pv.vz(),
    )
    return deltaR_vec_to_vec(sv_fd[1], gen_fd)


############################################################ Initialization ###
class PreWorker(fwliteworker.FwliteWorker):
    def node_setup(self, init_wrp):
        print "Starting:", init_wrp.sample, \
            init_wrp.event_handle.size(), init_wrp.filenames

        init_wrp.pre_worker = self
        self.init_wrp = init_wrp

        self.h_gen_particles = h_gen_particles
        self.h_pu_weight = h_pu_weight
        self.h_pv = h_pv
        self.h_ivf_normal = h_ivf_normal
        self.h_ivf_merged = h_ivf_merged
        self.h_jets = h_jets

        fs = self.result
        fs.NumInEvts = ROOT.TH1D('NumInEvts', 'NumInEvts', 1, 0.5, 1.5)
        fs.NumInEvts.Fill(1., init_wrp.event_handle.size())

    def node_process_event(self, event):
        self.weight = 1.
        event.getByLabel("offlinePrimaryVertices", self.h_pv)
        self.pv = self.h_pv.product()[0]
        event.getByLabel("selectedPatJetsPF", self.h_jets)
        self.jets = h_jets.product()

        # prepare ivf collections
        def mk_ivf_coll(label, handle):
            event.getByLabel(label, handle)
            coll = list(handle.product())
            for sv in coll:  # flight distances
                sv.__dict__.update(get_2d_3d_distances(self.pv, sv))
            coll = get_tuples_with_flight_dirs(coll, self.pv)
            setattr(self, label, coll)
        mk_ivf_coll("inclusiveMergedVerticesFiltered", self.h_ivf_normal)
        mk_ivf_coll("bToCharmDecayVertexMerged", self.h_ivf_merged)

        self.is_real_data = "Run" in self.init_wrp.sample
        if not self.is_real_data:
            event.getByLabel("genParticles", self.h_gen_particles)
            self.gen_particles = self.h_gen_particles.product()
            self.fin_b = final_b_hadrons(self.gen_particles)
            self.fin_d = final_d_hadrons(self.gen_particles)

            event.getByLabel("puWeight", "PUWeightTrue", self.h_pu_weight)
            self.weight *= self.h_pu_weight.product()[0]

    def node_finalize(self, init_wrp):
        del init_wrp.pre_worker
        print "Done:", init_wrp.sample, \
            init_wrp.event_handle.size(), init_wrp.filenames


############################################### Jet Plots for Normalization ###
class JetWorker(fwliteworker.FwliteWorker):
    def node_setup(self, init_wrp):
        self.init_wrp = init_wrp
        fs = self.result

        fs.make(
            "SelectedPatJetPt",
            ";selectedPatJets: p_{T}; number of events",
            100, 0., 500.
        )
        fs.make(
            "SelectedPatJetEta",
            ";selectedPatJets: #eta; number of events",
            100, -3.5, 3.5
        )

    def node_process_event(self, event):
        sample = self.init_wrp.sample
        if 'TTbar' in sample and 'NoMatch' not in sample:
            return

        fs = self.result
        pre_worker = self.init_wrp.pre_worker
        w = pre_worker.weight
        jets = pre_worker.jets
        for j in jets:
            fs.SelectedPatJetPt.Fill(j.pt(), w)
            fs.SelectedPatJetEta.Fill(j.eta(), w)


########################################### Vertex Selection and Histograms ###
class Worker(fwliteworker.FwliteWorker):
    def __init__(self, name, collection, filter_vtx=None):
        super(Worker, self).__init__(name)
        self.collection = collection
        self.n_matches_required = -1
        self.filter_vtx = filter_vtx

    def node_setup(self, init_wrp):
        self.init_wrp = init_wrp

        if 'TTbarNoMatch' == init_wrp.sample:
            self.n_matches_required = 0
        elif 'TTbarOneMatch' == init_wrp.sample:
            self.n_matches_required = 1
        elif 'TTbarTwoMatch' == init_wrp.sample:
            self.n_matches_required = 2
        elif 'TTbarBDMatch' == init_wrp.sample:
            self.n_matches_required = 3
        fs = self.result

        # Delta R plots
        fs.make(
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

        # Mass Template Fit
        fs.VertexMassVsDr = ROOT.TH2D(
            'VertexMassVsDr',
            ';Vertex #Delta R; Vertex Mass / GeV',
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

        # Mass Sum Template Fit
        fs.VertexMassSumVsDr = ROOT.TH2D(
            'VertexMassSumVsDr',
            ';Vertex #Delta R; Sum of Vertex Masses / GeV',
            100, 0., 5.,
            100, 0., 10.
        )
        fs.make(
            'VertexMassSumTemplate',
            ";B vertex candidate mass; number of events",
            100, 0., 10.
        )

        # Simultaneous Mass Template Fit
        fs.VertexMassPtLeadVsDr = ROOT.TH2D(
            'VertexMassPtLeadVsDr',
            ';Vertex #Delta R; Vertex Mass / GeV',
            100, 0., 5.,
            100, 0., 10.
        )
        fs.VertexMassPtSubLeadVsDr = ROOT.TH2D(
            'VertexMassPtSubLeadVsDr',
            ';Vertex #Delta R; Vertex Mass / GeV',
            100, 0., 5.,
            100, 0., 10.
        )

        # NTracks Template Fit
        fs.VertexNTracksVsDr = ROOT.TH2D(
            'VertexNTracksVsDr',
            ';Vertex #Delta R; Vertex number of tracks',
            100, 0., 5.,
            21, -.5, 20.5
        )
        fs.make(
            'VertexBeeNTracksTemplate',
            ";B vertex candidate number of tracks; number of events",
            21, -.5, 20.5
        )
        fs.make(
            'VertexDeeNTracksTemplate',
            ";D vertex candidate number of tracks; number of events",
            21, -.5, 20.5
        )

        # NTracks Sum Template Fit
        fs.VertexNTracksSumVsDr = ROOT.TH2D(
            'VertexNTracksSumVsDr',
            ';Vertex #Delta R; Vertex number of tracks',
            100, 0., 5.,
            21, -.5, 20.5
        )
        fs.make(
            'VertexNTracksSumTemplate',
            ";B vertex candidate number of tracks; number of events",
            21, -.5, 20.5
        )

        # Misc
        fs.make(
            'VtxNSharedTracks',
            ';number of shared tracks; number of events',
            11, -.5, 10.5
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
            ";pt leading vertex: no. tracks (weight > 0.5); number of events",
            21, -.5, 20.5
        )
        fs.make(
            "VtxPtSubLeadNumTracks",
            ";pt sub leading vertex: no. tracks (weight > 0.5); number of events",
            21, -.5, 20.5
        )
        fs.make(
            "VtxPtLeadTracksSize",
            ";pt leading vertex: no. tracks; number of events",
            21, -.5, 20.5
        )
        fs.make(
            "VtxPtSubLeadTracksSize",
            ";pt sub leading vertex: no. tracks; number of events",
            21, -.5, 20.5
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
            ";Vertex Mass / GeV;number of vertices",
            100, 0., 10.
        )
        fs.make(
            "VtxDeeMass",
            ";Vertex Mass / GeV;number of vertices",
            100, 0., 10.
        )
        fs.make(
            "VtxBeeDeeMatchSig",
            ";#DeltaX for matching B and D;number of vertices",
            100, 0., 100.
        )
        fs.make(
            "VertexBeeDistLtDeeDist",
            ";Bee.dist3d() < Dee.dist3d() ;number of events",
            2, -0.5, 1.5
        )
        fs.make(
            "VertexBeeMassLtDeeMass",
            ";Bee.p4().mass() < Dee.p4().mass() ;number of events",
            2, -0.5, 1.5
        )
        fs.make(
            "VertexBeeDistLtDeeDistOneSigma",
            ";Bee.dist3d() < Dee.dist3d() ;number of events",
            2, -0.5, 1.5
        )
        fs.VertexBeeVsDeeDist2D = ROOT.TH2D(
            'VertexBeeVsDeeDist2D',
            ';D Vertex Dist2D; B Vertex Dist2D',
            100, 0., 10.,
            100, 0., 10.
        )
        fs.VertexBeeVsDeeDist3D = ROOT.TH2D(
            'VertexBeeVsDeeDist3D',
            ';D Vertex Dist3D; B Vertex Dist3D',
            100, 0., 10.,
            100, 0., 10.
        )
        fs.VertexBeeVsDeeSig2D = ROOT.TH2D(
            'VertexBeeVsDeeSig2D',
            ';D Vertex Sig2D; B VertexSig2D',
            100, 0., 10.,
            100, 0., 10.
        )
        fs.VertexBeeVsDeeSig3D = ROOT.TH2D(
            'VertexBeeVsDeeSig3D',
            ';D Vertex Sig3D; B Vertex Sig3D',
            100, 0., 10.,
            100, 0., 10.
        )

    def node_process_event(self, event):
        fs = self.result
        pre_worker = self.init_wrp.pre_worker

        # ivf vertices
        ivf_vtx_fd = getattr(pre_worker, self.collection)
        if self.filter_vtx:
            ivf_vtx_fd = self.filter_vtx(ivf_vtx_fd)
        ivf_vtx_fd_pt_sort = sorted(ivf_vtx_fd, key=lambda v: -v[0].p4().pt())

        # if more than two ivf vtx, take the ones closest together
        ivf_vtx_fd_max2 = ivf_vtx_fd
        if len(ivf_vtx_fd) > 2:
            ivf_vtx_fd_max2 = min(
                itertools.combinations(ivf_vtx_fd, 2),
                key=lambda c: fd_dr(*c)
            )

        # matching
        matched = []
        matched_bd = []
        fin_b = []
        fin_bd = []
        if not pre_worker.is_real_data:
            gen_particles = pre_worker.gen_particles
            fin_b = pre_worker.fin_b

            matched = matching(
                ivf_vtx_fd_max2,
                fin_b,
                #lambda a, b: deltaR_vec_to_cand(a[1], b),
                #lambda a, b: deltaR_cand_to_cand(a[0], b),
                lambda a, b: deltaR_from_pv(a, b, fin_b[0]),
                DR_for_matching
            )

            if len(matched) == 1:
                # try to match D hadron as well
                fin_b_vtx, fin_b_gp, _ = matched[0]
                fin_d = pre_worker.fin_d  # final_d_hadrons(get_all_daughters(gen_particles, [fin_b_gp]))
                other_vtxs = list(v for v in ivf_vtx_fd_max2 if v != fin_b_vtx)
                if other_vtxs:
                    matched_d = matching(
                        other_vtxs,
                        fin_d,
                        #lambda a, b: deltaR_cand_to_cand(a[0], b),
                        lambda a, b: deltaR_from_pv(a, b, fin_b[0]),
                        DR_for_matching
                    )
                    matched_bd = matched + matched_d

        n_matched = len(matched)
        if 1 == n_matched and 2 == len(matched_bd):
            n_matched = 3

        ################################################### fill histograms ###

        # check for right ttbar or different sample
        if self.n_matches_required not in (-1, n_matched):
            return

        w = pre_worker.weight
        fs.NumIvfVertices.Fill(len(ivf_vtx_fd), w)

        if ivf_vtx_fd_pt_sort:
            fs.VtxPtLeadMass.Fill(ivf_vtx_fd_pt_sort[0][0].p4().mass(), w)
            fs.VtxPtLeadPt.Fill(ivf_vtx_fd_pt_sort[0][0].p4().pt(), w)
            fs.VtxPtLeadEta.Fill(ivf_vtx_fd_pt_sort[0][0].p4().eta(), w)
            fs.VtxPtLeadFdEta.Fill(ivf_vtx_fd_pt_sort[0][1].eta(), w)
            fs.VtxPtLeadNumTracks.Fill(ivf_vtx_fd_pt_sort[0][0].nTracks(), w)
            fs.VtxPtLeadTracksSize.Fill(
                ivf_vtx_fd_pt_sort[0][0].tracksSize(), w)
        if len(ivf_vtx_fd_pt_sort) > 1:
            fs.VtxPtSubLeadMass.Fill(ivf_vtx_fd_pt_sort[1][0].p4().mass(), w)
            fs.VtxPtSubLeadPt.Fill(ivf_vtx_fd_pt_sort[1][0].p4().pt(), w)
            fs.VtxPtSubLeadEta.Fill(ivf_vtx_fd_pt_sort[1][0].p4().eta(), w)
            fs.VtxPtSubLeadFdEta.Fill(ivf_vtx_fd_pt_sort[1][1].eta(), w)
            fs.VtxPtSubLeadNumTracks.Fill(ivf_vtx_fd_pt_sort[1][0].nTracks(), w)
            fs.VtxPtSubLeadTracksSize.Fill(
                ivf_vtx_fd_pt_sort[1][0].tracksSize(), w)

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
            if hasattr(ivf_vtx_fd_max2[0][0], 'originalMass'):
                mass0 = ivf_vtx_fd_max2[0][0].originalMass()
                mass1 = ivf_vtx_fd_max2[1][0].originalMass()
                n_trk0 = ivf_vtx_fd_max2[0][0].originalNTracks()
                n_trk1 = ivf_vtx_fd_max2[1][0].originalNTracks()
                mass_lead = ivf_vtx_fd_pt_sort[0][0].originalMass()
                mass_subl = ivf_vtx_fd_pt_sort[1][0].originalMass()
            else:
                mass0 = ivf_vtx_fd_max2[0][0].p4().mass()
                mass1 = ivf_vtx_fd_max2[1][0].p4().mass()
                n_trk0 = ivf_vtx_fd_max2[0][0].nTracks()
                n_trk1 = ivf_vtx_fd_max2[1][0].nTracks()
                mass_lead = ivf_vtx_fd_pt_sort[0][0].p4().mass()
                mass_subl = ivf_vtx_fd_pt_sort[1][0].p4().mass()
            fs.VertexMassVsDr.Fill(dr_mom, mass0, w)
            fs.VertexMassVsDr.Fill(dr_mom, mass1, w)
            fs.VertexMassSumVsDr.Fill(dr_mom, mass0 + mass1, w)
            fs.VertexNTracksVsDr.Fill(dr_mom, n_trk0, w)
            fs.VertexNTracksVsDr.Fill(dr_mom, n_trk1, w)
            fs.VertexNTracksSumVsDr.Fill(dr_mom, n_trk0 + n_trk1, w)
            fs.VertexMassPtLeadVsDr.Fill(dr_mom, mass_lead, w)
            fs.VertexMassPtSubLeadVsDr.Fill(dr_mom, mass_subl, w)

            # fit templates
            if dr_mom < 0.1:
                #b_cand, d_cand = sorted(ivf_vtx_fd_max2, key=lambda v: v[0].dxyz_val)
                #if (b_cand[0].dxyz_val + b_cand[0].dxyz_err <
                #        d_cand[0].dxyz_val - d_cand[0].dxyz_err):
                b_cand, d_cand = sorted(ivf_vtx_fd_max2, key=lambda v: -v[0].p4().mass())
                if True:
                    b_mass = b_cand[0].p4().mass()
                    b_ntrk = b_cand[0].nTracks()
                    d_mass = d_cand[0].p4().mass()
                    d_ntrk = d_cand[0].nTracks()
                    fs.VertexBeeMassTemplate.Fill(b_mass, w)
                    fs.VertexDeeMassTemplate.Fill(d_mass, w)
                    fs.VertexMassSumTemplate.Fill(b_mass + d_mass, w)
                    fs.VertexBeeNTracksTemplate.Fill(b_ntrk, w)
                    fs.VertexDeeNTracksTemplate.Fill(d_ntrk, w)
                    fs.VertexNTracksSumTemplate.Fill(b_ntrk + d_ntrk, w)


            # set comparator for tracks
            ROOT.reco.Track.__eq__ = lambda a, b: \
                round(a.pt(), 10) == round(b.pt(), 10) \
                and round(a.eta(), 10) == round(b.eta(), 10)
            ROOT.reco.Track.__hash__ = lambda a: \
                hash(round(a.pt() + a.eta(), 10))

            # n shared tracks
            all_tracks = list(t
                              for v in ivf_vtx_fd_max2
                              for t in v[0].refittedTracks())
            unique_tracks = set(all_tracks)
            fs.VtxNSharedTracks.Fill(len(all_tracks) - len(unique_tracks), w)

        ######################################################### MC histos ###
        if not pre_worker.is_real_data:
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
            bee, bee_gen, bee_match_dr = bee
            dee, dee_gen, dee_match_dr = dee
            bee, dee = bee[0], dee[0]
            fs.VtxBeeNumTracks.Fill(bee.nTracks(), w)
            fs.VtxDeeNumTracks.Fill(dee.nTracks(), w)
            fs.VtxBeeMass.Fill(bee.p4().M(), w)
            fs.VtxDeeMass.Fill(dee.p4().M(), w)
            vtx_dr = deltaR_cand_to_cand(bee, dee)
            if bee_match_dr < vtx_dr < 0.1 and dee_match_dr < vtx_dr:
                fs.VertexBeeVsDeeDist3D.Fill(dee.dxyz_val, bee.dxyz_val, w)
                fs.VertexBeeVsDeeDist2D.Fill(dee.dxy_val, bee.dxy_val, w)
                fs.VertexBeeVsDeeSig3D.Fill(dee.dxyz_val/dee.dxyz_err, bee.dxyz_val/bee.dxyz_err, w)
                fs.VertexBeeVsDeeSig2D.Fill(dee.dxy_val/dee.dxy_err, bee.dxy_val/bee.dxy_err, w)
                fs.VertexBeeDistLtDeeDist.Fill(dee.dxyz_val > bee.dxyz_val, w)
                fs.VertexBeeMassLtDeeMass.Fill(dee.p4().M() > bee.p4().M, w)
                if (dee.dxyz_val-dee.dxyz_err) > (bee.dxyz_val+bee.dxyz_err):
                    fs.VertexBeeDistLtDeeDistOneSigma.Fill(dee.dxyz_val > bee.dxyz_val, w)

            # matching significance
            matching_bd_cov = matching(
                ivf_vtx_fd_max2,
                fin_bd,
                lambda a, b: covariance_significance(a[0], b),
                100.
            )
            for _, _, sig in matching_bd_cov:
                fs.VtxBeeDeeMatchSig.Fill(sig, w)


############################################################ Vertex Filters ###
vtx_cuts = lambda vtx: filter(lambda v: (
    abs(v[0].p4().eta()) < 2
    and v[0].p4().pt() > 8
#    and 6.5 > v[0].p4().mass() > 1.4
    and v[0].nTracks() > 2
    and deltaR_vec_to_vec(v[1], v[0].p4()) < 0.1
    and v[0].dxy_val < 2.5
    and (v[0].dxy_val / v[0].dxy_err) > 3.
    and (v[0].dxyz_val / v[0].dxyz_err) > 5.
), vtx)


has_dr_lt = lambda vtx_fd, bound: any(
    deltaR_cand_to_cand(a[0], b[0]) < bound
    for a, b in itertools.combinations(vtx_fd, 2)
)
vtx_dr_lt_0p2 = lambda vtx: vtx if has_dr_lt(vtx, 0.2) else list()
vtx_no_dr_lt_1p0 = lambda vtx: vtx if not has_dr_lt(vtx, 1.0) else list()


####################################################### Worker Organization ###
workers = [
    PreWorker("PreWorker"),
    JetWorker("JetWorker"),
    #Worker("IvfMerged", "inclusiveMergedVertices"),
    Worker("IvfMergedFilt", "inclusiveMergedVerticesFiltered"),
    #Worker("IvfMergedFiltLt0p2", 'inclusiveMergedVerticesFiltered', vtx_dr_lt_0p2),
    #Worker("IvfMergedFiltGt1p0", 'inclusiveMergedVerticesFiltered', vtx_no_dr_lt_1p0),
    Worker("IvfMergedFiltCuts", "inclusiveMergedVerticesFiltered", vtx_cuts),
    Worker("IvfB2cMerged", "bToCharmDecayVertexMerged"),
    Worker("IvfB2cMergedCuts", "bToCharmDecayVertexMerged", vtx_cuts),
]


if __name__ == '__main__':
    fwliteworker.work(workers)


