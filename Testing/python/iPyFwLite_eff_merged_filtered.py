# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <markdowncell>

# # IVF BTagDR Efficiency
# ### imports, parameters, functions

# <codecell>

import ROOT
import itertools
import glob

from rootpy.plotting import Hist
from MyUtility.PythonUtil.genParticles import final_b_mesons
from DataFormats.FWLite import Events,Handle
#from dcache_TTdilep import files

# <codecell>

# parameters
DR_for_matching = 0.1
files = map(lambda f: "file:%s"%f, glob.glob(
    "/nfs/dust/cms/user/tholenhe/samples/"
    "DoubleVtxEffTTdilep/TTdilep_presel_*.root"
))

# <codecell>

# cmssw 
h_genParticles = Handle("vector<reco::GenParticle>")
h_pv = Handle("vector<reco::Vertex>")
h_ivf = Handle("vector<reco::Vertex>")
events = iter(Events(files))

# functions
DeltaR = ROOT.Math.VectorUtil.DeltaR
deltaR = lambda a,b: DeltaR(a.p4(), b.p4())
def my_deltaR(a, b):
    return (
        (a.eta()-b.eta())**2 
        + (float(a.phi())-float(b.phi()))**2
    )**.5

def get_all_flight_dirs(vertices, primary_vertex):
    def mkfd(vtx):
        pv, sv = primary_vertex, vtx.position()
        return ROOT.Math.XYZVector(
            sv.x() - pv.x(),
            sv.y() - pv.y(),
            sv.z() - pv.z(),
        )
    return list(mkfd(sv) for sv in vertices)

def get_matched_flightdirs(flightdirs, gen_final_bs):
    unmatched_final_bs = gen_final_bs[:]
    def is_close_to_fin_b(fd):
        for bee in unmatched_final_bs[:]:
            if my_deltaR(fd, bee.p4()) < DR_for_matching:
                unmatched_final_bs.remove(bee)  ## not to be matched twice
                return True
        return False            
    
    return list(
        fd for fd in flightdirs
        if is_close_to_fin_b(fd)
    )

# <markdowncell>

# ### histograms

# <codecell>

from varial import diskio
fs = diskio.fileservice("fs_ivf_merged_filtered")

# all events
fs.NumFinalBs = ROOT.TH1D(
    "NumFinalBs",
    ";number of final B's;number of events",
    8, -.5, 7.5
)
fs.NumIvfVertices = ROOT.TH1D(
    "NumIvfVertices",
    ";number of IVF vertices;number of events",
    8, -.5, 7.5
)

# with n_ivf == 2, no / one / two matched
fs.VertexDR = ROOT.TH1D(
    "VertexDR",
    ";#Delta R;number of vertices",
    100, 0., 5.
)
fs.VertexDRNoMatch = ROOT.TH1D(
    "VertexDRNoMatch",
    ";#Delta R;number of vertices",
    100, 0., 5.
)
fs.VertexDROneMatch = ROOT.TH1D(
    "VertexDROneMatch",
    ";#Delta R;number of vertices",
    100, 0., 5.
)
fs.VertexDRTwoMatch = ROOT.TH1D(
    "VertexDRTwoMatch",
    ";#Delta R;number of vertices",
    100, 0., 5.
)

# <markdowncell>

# ### eventloop

# <codecell>

vtx_dr_histos = [fs.VertexDRNoMatch, fs.VertexDROneMatch, fs.VertexDRTwoMatch]

print "starting eventloop"
for event in events:

    # final B hadron generator particles
    event.getByLabel("genParticles", h_genParticles)
    genParticles = h_genParticles.product()
    fin_bs = final_b_mesons(genParticles)

    # ivf vertices

    event.getByLabel("goodOfflinePrimaryVertices", h_pv)
    event.getByLabel("inclusiveMergedVerticesFiltered", h_ivf)
    flightdirs = get_all_flight_dirs(h_ivf.product(), h_pv.product()[0])
    matched_fds = get_matched_flightdirs(flightdirs, fin_bs)
    
    # fill histos for all events
    fs.NumFinalBs.Fill(len(fin_bs))
    fs.NumIvfVertices.Fill(len(flightdirs))
    
    # fill histos for n_matched == 2
    if len(flightdirs) == 2:
        flightdir_dR = my_deltaR(*flightdirs)
        fs.VertexDR.Fill(flightdir_dR)
        vtx_dr_histos[len(matched_fds)].Fill(flightdir_dR)


print "eventloop done"

