import itertools
import os

import ROOT
ROOT.gROOT.SetBatch()
import varial
import varial.main
import varial.sample


tmp_dataset = '_tmp_das_datasets.txt'
tmp_files = '_tmp_das_files_%s.txt'
dbs_command_datasets = 'das_client.py --limit 10000 --query "' \
    'instance=prod/phys03 ' \
    'dataset=/*/*tholen-TTdil*/USER' \
    '" | grep htholen > ' + tmp_dataset
dbs_command_files = 'das_client.py --limit 10000 --query "' \
    'instance=prod/phys03 ' \
    'file dataset=%s' \
    '" | grep htholen > ' + tmp_files
output_dir = '/nfs/dust/cms/user/tholenhe/samples/TTdilep/'


def parse_sample_name(dataset):
    #e.g. /MuEG/htholen-TTdilepData_RunA-d4aa20010c5c4332ad27156a5b618077/USER
    #                               ^^^^
    if 'TTdilepAllMC' in dataset:
        back = dataset.split('TTdilepAllMC_')[1]
    else:
        back = dataset.split('TTdilepData_')[1]
    return back.split('-')[0]


def main():
    used_das_client = False

    # query datasets
    if not os.path.exists(tmp_dataset):
        os.system(dbs_command_datasets)
        used_das_client = True
    with open('_tmp_das_datasets.txt') as f:
        datasets = f.readlines()
    datasets = dict((parse_sample_name(d), d) for d in datasets)

    # query files
    files = {}
    for sample, dataset in datasets.iteritems():
        if not os.path.exists(tmp_files % sample):
            os.system(dbs_command_files % (dataset, sample))
            used_das_client = True
        with open(tmp_files % sample) as f:
            files[sample] = f.readlines()

    if used_das_client:
        os.system("head -3 _tmp_das_*")
        print "WARNING I used the das client to query for filenames."
        print "WARNING Therefore, I do not execute cmsRun. Please make "
        print "WARNING sure the _tmp_das_* files (see above) are not empty."
        return

    samples = {}
    for sample, files in files.iteritems():

        # group with number of files
        if "TTbar" == sample:
            groups = itertools.izip(*[iter(files)]*5)
        else:
            groups = itertools.izip(*[iter(files)]*20)

        # make samples
        for i, grp in enumerate(groups):
            name = '%s_%d' % (sample, i)
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
        suppress_cmsRun_exec=used_das_client,
        cfg_main_import_path="BTagDeltaR.PreSelections.MergeTtdilep_cfg"
    )


if __name__ == '__main__':
    main()