import itertools
import os

import ROOT
ROOT.gROOT.SetBatch()
import varial
import varial.main
import varial.sample


dbs_command_datasets = 'das_client.py --query "' \
    'instance=prod/phys03 ' \
    'dataset=/*/htholen-TTdil*/USER' \
    '" | grep htholen > _tmp_das_datasets.txt'
dbs_command_files = 'das_client.py --query "' \
    'instance=prod/phys03 ' \
    'file dataset=%s' \
    '" --limit 10000 | grep htholen > _tmp_das_files_%s.txt'
output_dir = '/nfs/dust/cms/user/tholenhe/samples/TTdilep/'


def parse_sample_name(dataset):
    #e.g. /MuEG/htholen-TTdilepData_RunA-d4aa20010c5c4332ad27156a5b618077/USER
    #                               ^^^^
    if 'TTdilepAllMC' in dataset:
        back = dataset.split['TTdilepAllMC_'][1]
    else:
        back = dataset.split['TTdilepData_'][1]
    return back.split['-'][0]


def main():
    # query datasets
    os.system(dbs_command_datasets)
    with open('_tmp_das_datasets.txt') as f:
        datasets = f.readlines()
    datasets = dict((parse_sample_name(d), d) for d in datasets)

    # query files
    files = {}
    for sample, dataset in datasets.iteritems():
        os.system(dbs_command_files % (dataset, sample))
        with open('_tmp_das_files_%s.txt' % sample) as f:
            files[sample] = f.readlines()

    samples = {}
    for sample, files in files.iteritems():

        # group with number of files
        groups = itertools.izip(*[iter(files)]*10)

        # make samples
        for i, grp in enumerate(groups):
            name = '%s_%i' % sample, i
            samples[name] = varial.sample.Sample(
                name=name,
                is_data="TTdilepData" in grp[0],
                input_files=grp,
                output_file=output_dir,
                lumi=1.
            )

    varial.main.main(
        samples=samples,
        cfg_use_file_service=False,
        try_reuse_results=True,
        suppress_cmsRun_exec=True,
        cfg_main_import_path="BTagDeltaR.PreSelections.MergeTtdilep_cfg"
    )