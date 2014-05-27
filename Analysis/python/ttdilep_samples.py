from varial.sample import Sample
from glob import glob


class TTDilepEMuSample(Sample):
    def __init__(self, **kws):
        kws['input_files'] = glob(
            '/nfs/dust/cms/user/tholenhe/samples/TTdilep/'
            '%s_*.root' % kws['name']
        )

        # with open('_tmp_das_files_%s.txt' % kws['name']) as f:
        #     kws['input_files'] = map(lambda s:
        #         'dcap://dcache-cms-dcap.desy.de/pnfs/desy.de/cms/tier2'
        #         + s.strip(),
        #         f.readlines()
        #     )

        super(TTDilepEMuSample, self).__init__(**kws)


smp_emu = list()


######################################################################## MC ###
smp_emu.append(TTDilepEMuSample(
    name='DiBosonWW',
    x_sec=57.11,
    n_events=10000431,
    legend='Di-boson'
))

smp_emu.append(TTDilepEMuSample(
    name='DiBosonWZ',
    x_sec=32.32,
    n_events=10000283,
    legend='Di-boson',
))

smp_emu.append(TTDilepEMuSample(
    name='DiBosonZZ',
    x_sec=8.26,
    n_events=9752010,
    legend='Di-boson',
))

smp_emu.append(TTDilepEMuSample(
    name='SingleTW',
    x_sec=11.1,
    n_events=497658,
    legend='Single top',
))

smp_emu.append(TTDilepEMuSample(
    name='SingleTbarW',
    x_sec=11.1,
    n_events=493460,
    legend='Single top',
))

smp_emu.append(TTDilepEMuSample(
    name='ZJets10to50',
    x_sec=11050.0,
    n_events=37835275,
    legend='DY + Jets',
))

smp_emu.append(TTDilepEMuSample(
    name='ZJets50plus',
    x_sec=2950.0,
    n_events=30458871,
    legend='DY + Jets',
))

smp_emu.append(TTDilepEMuSample(
    name='WJets',
    x_sec=37509.,
    n_events=57709905,
    legend='DY + Jets',
))

smp_emu.append(TTDilepEMuSample(
    name='TTbar',
    x_sec=245.,
    n_events=12116717,
    legend='TTbar',
))


###################################################################### Data ###
# smp_emu.append(TTDilepEMuSample(
#     name='RunA',
#     is_data=True,
#     lumi=94.676 + 390.113 + 391.376,  # pixelLumiCalc.py
#     legend='Data',
# ))
#
# smp_emu.append(TTDilepEMuSample(
#     name='RunB',
#     is_data=True,
#     lumi=4411.0,  # pixelLumiCalc.py
#     legend='Data',
# ))
#
# smp_emu.append(TTDilepEMuSample(
#     name='RunC',
#     is_data=True,
#     lumi=1734. + 5321.,  # pixelLumiCalc.py
#     legend='Data',
# ))
#
# smp_emu.append(TTDilepEMuSample(
#     name='RunD',
#     is_data=True,
#     lumi=7360.,  # pixelLumiCalc.py
#     legend='Data',
# ))

