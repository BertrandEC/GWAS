import pandas as pd
import tempfile
import subprocess
import os
from enum import Enum


class LDPopulation(Enum):
    EUROPEAN = 'https://zenodo.org/records/10515792/files/1000G_Phase3_baselineLD_v2.2_ldscores.tgz?download=1'
    EAST_ASIAN = 'https://zenodo.org/records/10515792/files/1000G_Phase3_EAS_baselineLD_v2.2_ldscores.tgz?download=1'


def ldsc(disease1: str, n_disease1, disease2: str, n_disease2):
    with tempfile.TemporaryDirectory(delete=False) as tmpdir:
        print(f'{tmpdir}\n\n')

        munge_cmd1 = f'munge_sumstats.py --sumstats {disease1} --N {n_disease1} --snp variant_id --ignore rsid --out {tmpdir}/disease1' # --merge-alleles data/hm3_no_MHC.list.txt'
        munge_cmd2 = f'munge_sumstats.py --sumstats {disease2} --N {n_disease2} --snp variant_id --ignore rsid --out {tmpdir}/disease2' # --merge-alleles data/hm3_no_MHC.list.txt'
        subprocess.run(munge_cmd1, shell=True)
        subprocess.run(munge_cmd2, shell=True)
        cmd = f'ldsc.py --rg {tmpdir}/disease1.sumstats.gz,{tmpdir}/disease2.sumstats.gz --ref-ld-chr data/eur_w_ld_chr/ --w-ld-chr data/eur_w_ld_chr/ --out {tmpdir}/ldsc_results'
        print(os.getcwd())
        subprocess.run(cmd, shell=True)
        print(f'{tmpdir}/ldsc_results')