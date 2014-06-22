import ROOT

from varial.sample import Sample
from varial import settings
from glob import glob


settings.colors = {
    'Di-boson': ROOT.kYellow,
    'Single top': ROOT.kMagenta,
    'DY + jets': ROOT.kBlue,
    'W + jets': ROOT.kOrange,
    'TTbar': ROOT.kRed,
    'TTbar (no match)': ROOT.kRed + 3,
    'TTbar (1 match)': ROOT.kRed + 1,
    'TTbar (B+D match)': ROOT.kRed - 4,
    'TTbar (2 matches)': ROOT.kRed - 9,
}
settings.stacking_order = [
    'Di-boson',
    'Single top',
    'DY + jets',
    'W + jets',
    'TTbar',
    'TTbar (no match)',
    'TTbar (1 match)',
    'TTbar (B+D match)',
    'TTbar (2 matches)',
]


class TTDilepEMuSample(Sample):
    def __init__(self, **kws):
        filename = kws['name']
        if 'TTbar' in filename:
            filename = 'TTbar'
        kws['input_files'] = glob(
            '/nfs/dust/cms/user/tholenhe/samples/TTdilep/'
            '%s_*.root' % filename
        )

        # with open('_tmp_das_files_%s.txt' % kws['name']) as f:
        #     kws['input_files'] = map(lambda s:
        #         'dcap://dcache-cms-dcap.desy.de/pnfs/desy.de/cms/tier2'
        #         + s.strip(),
        #         f.readlines()
        #     )

        super(TTDilepEMuSample, self).__init__(**kws)


######################################################################## MC ###
smp_emu_mc = list()


smp_emu_mc.append(TTDilepEMuSample(
    name='DiBosonWW',
    x_sec=57.11,
    n_events=10000431,
    legend='Di-boson'
))

smp_emu_mc.append(TTDilepEMuSample(
    name='DiBosonWZ',
    x_sec=32.32,
    n_events=10000283,
    legend='Di-boson',
))

smp_emu_mc.append(TTDilepEMuSample(
    name='DiBosonZZ',
    x_sec=8.26,
    n_events=9752010,
    legend='Di-boson',
))

smp_emu_mc.append(TTDilepEMuSample(
    name='SingleTW',
    x_sec=11.1,
    n_events=497658,
    legend='Single top',
))

smp_emu_mc.append(TTDilepEMuSample(
    name='SingleTbarW',
    x_sec=11.1,
    n_events=493460,
    legend='Single top',
))

smp_emu_mc.append(TTDilepEMuSample(
    name='ZJets10to50',
    x_sec=11050.0,
    n_events=37835275,
    legend='DY + jets',
))

smp_emu_mc.append(TTDilepEMuSample(
    name='ZJets50plus',
    x_sec=2950.0,
    n_events=30458871,
    legend='DY + jets',
))

smp_emu_mc.append(TTDilepEMuSample(
    name='WJets',
    x_sec=37509.,
    n_events=57709905,
    legend='W + jets',
))

smp_emu_mc.append(TTDilepEMuSample(
    name='TTbarNoMatch',
    x_sec=245./9.,
    n_events=12116717,  # * 15051. / 522465.,
    legend='TTbar (no match)',
))

smp_emu_mc.append(TTDilepEMuSample(
    name='TTbarOneMatch',
    x_sec=245./9.,
    n_events=12116717,  # * 15051. / 522465.,
    legend='TTbar (1 match)',
))

smp_emu_mc.append(TTDilepEMuSample(
    name='TTbarBDMatch',
    x_sec=245./9.,
    n_events=12116717,  # * 15051. / 522465.,
    legend='TTbar (B+D match)',
))

smp_emu_mc.append(TTDilepEMuSample(
    name='TTbarTwoMatch',
    x_sec=245./9.,
    n_events=12116717,  # * 15051. / 522465.,
    legend='TTbar (2 matches)',
))

"""
smp_emu_mc.append(TTDilepEMuSample(
    name='TTbar',
    x_sec=245./9.,
    n_events=12116717,  # * 15051. / 522465.,
    legend='TTbar',
))
"""
# presel: 522465
# smp_emu_mc[-1].input_files = smp_emu_mc[-1].input_files[:12]

###################################################################### Data ###
smp_emu_data = list()


smp_emu_data.append(TTDilepEMuSample(
    name='RunA',
    is_data=True,
    lumi=94.676 + 390.113 + 391.376,  # pixelLumiCalc.py
    legend='Data',
))

smp_emu_data.append(TTDilepEMuSample(
    name='RunB',
    is_data=True,
    lumi=4411.0,  # pixelLumiCalc.py
    legend='Data',
))

smp_emu_data.append(TTDilepEMuSample(
    name='RunC',
    is_data=True,
    lumi=1734. + 5321.,  # pixelLumiCalc.py
    legend='Data',
))

smp_emu_data.append(TTDilepEMuSample(
    name='RunD',
    is_data=True,
    lumi=7360.,  # pixelLumiCalc.py
    legend='Data',
))

