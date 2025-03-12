import tempfile
import subprocess
import re
from enum import Enum
from dataclasses import dataclass
from pathlib import Path

class LDSCError(Exception):
    pass

class LDPopulation(Enum):
    EUROPEAN = ('eur_w_ld_chr', 'https://ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/')
    EAST_ASIAN = (None, None)

ALLELE_LIST = ('hm3_no_MHC.list.txt', 'https://ibg.colorado.edu/cdrom2021/Day06-nivard/GenomicSEM_practical/')

@dataclass
class LDSCResult:
     correlation: float
     correlation_std_error: float
     p_value: float


def ldsc(disease1: str, n_disease1: int, disease2: str, n_disease2: int, population: LDPopulation = LDPopulation.EUROPEAN) -> LDSCResult | None:
    """Returns results from running LDSC on given"""
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            if not Path("data", ALLELE_LIST[0]).is_file():
                raise LDSCError(f'Allele list file not found at data/{ALLELE_LIST[0]}. Please download from {ALLELE_LIST[1]}')
            munge_cmd = f'munge_sumstats.py --sumstats {{}} --N {{}} --snp variant_id --ignore rsid --out {tmpdir}/disease1 --merge-alleles data/{ALLELE_LIST[0]}'
            subprocess.run(munge_cmd.format(disease1, n_disease1), shell=True, capture_output=True, check=True)
            subprocess.run(munge_cmd.format(disease2, n_disease2), shell=True, capture_output=True, check=True)
            if not Path("data", population.value[0]).exists():
                raise LDSCError(f'Precomputed LDScores do not exist at data/{population.value[0]}. Please download from {population.value[1]}')
            cmd = f'ldsc.py --rg {tmpdir}/disease1.sumstats.gz,{tmpdir}/disease2.sumstats.gz --ref-ld-chr data/{population.value[0]}/ --w-ld-chr data/{population.value[0]}/ --out {tmpdir}/ldsc_results'
            full_result = subprocess.run(cmd, shell=True, capture_output=True, check=True).stdout.decode('utf-8')
            ldsc_results =  full_result.partition("Genetic Correlation")[2]
            match = re.search(r"^Genetic Correlation: (.+) \((.+)\)\s.*\s^P: (.+)", ldsc_results, re.MULTILINE)
            if not match: return None
            (correlation, std_error, p_value) = match.groups()
            return LDSCResult(correlation, std_error, p_value)
        except subprocess.CalledProcessError:
          return None