# conda 환경: conda activate cs2

import sys, os
import subprocess
import pandas as pd
from glob import glob
import shutil


class ABL1VUS:
    def __init__(self, sample_id:str, r1, r2, exon:str):
        '''CML VUS screening에서 ABL1의 variants counting을 위한 함수
        CRISPResso를 이용한 counting을 하고, 그 output을 이용해 variants list별로 정리된 
        output 파일을 생성한다. '''
        
        self._check_input(exon)

        self.info = self._exon_align_info(exon)
        
        self.sample_id = sample_id.replace('.', '_')
        self.exon = exon

        self.refseq = self.info['refseq']
        self.center = self.info['center']
        self.window = self.info['window']

        data    = f'-r1 {r1} -r2 {r2} -n {self.sample_id}'
        align   = f'-a {self.refseq} -an {exon} -g {self.center} --plot_window_size {self.window}'
        output  = f'--file_prefix {self.sample_id}'
        
        self.command = f'CRISPResso {data} {align} {output}'


    def run(self, out_dir:str, save_plot:bool=False, remove_temp=True):

        command  = self.command + f' -o {out_dir}'

        if save_plot == False:
            command = command + f' --suppress_plots --suppress_report'

        # Run CRISPResso
        subprocess.run([command], shell=True, check=True)

        # Copy and rename frequency table
        file_from = f'{out_dir}/CRISPResso_on_{self.sample_id}/{self.sample_id}.{self.exon}.Alleles_frequency_table_around_sgRNA_{self.center}.txt'
        file_to   = f'{out_dir}/NGS_frequency_table/{self.sample_id}.txt'

        shutil.copyfile(file_from, file_to) 

        if remove_temp == True:
            self._remove_temp_files(out_dir)

    def _remove_temp_files(self, out_dir:str):
        shutil.rmtree(f'{out_dir}/CRISPResso_on_{self.sample_id}')



    def _check_input(self, exon):
        if exon not in ['exon4', 'exon5', 'exon6', 'exon7', 'exon8', 'exon9', 'invivo_exon4', 'invivo_exon9']:
            raise ValueError('Not available exon input. Exon list: exon4, exon5, exon6, exon7, exon8, exon9')
    

    def _exon_align_info(self, exon):
        # Refseq을 인식해서 center sequence를 자동으로 만들고, 그 주위로 window 길이도 정해주게 만들 수 있을듯. 
        # 코드 구현이 되면 ABL1으로만 제한되지 않고, refseq을 input으로 받아서 모든 종류의 gene에 대해 사용할 수 있는 함수가 된다. 
        # Reference 정보 정리하면서 코드 수정해보기 

        dict_abl1_info = {

            'exon4': {
                'refseq': 'ttgagcttgcctgtctctgtgggctgaaggctgttccctgtttccttcagctctacgtctcctccgagagccgcttcaacaccctggccgagttggttcatcatcattcaacggtggccgacgggctcatcaccacgctccattatccagccccaaagcgcaacaagcccactgtctatggtgtgtcccccaactacgacaagtgggagatggaacgcacggacatcaccatgaagcacaagctgggcgggggccagtacggggaggtgtacgagggcgtgtggaagaaatacagcctgacggtggccgtgaagaccttgaaggtaggctgggactgccgggggtgcccagggtacgtggggcaag',
                'center': 'aagcccactgtctatggtgt',
                'window': 163,
            },
            'exon4_150PE_1': {
                'refseq': 'ttgagcttgcctgtctctgtgggctgaaggctgttccctgtttccttcagctctacgtctcctccgagagccgcttcaacaccctggccgagttggttcatcatcattcaacggtggccgacgggctcatcaccacgctccattatccagccccaaagcgcaacaagcccactgtctatggtgtgtcccccaactacgacaagtgggagatggaacgcacggac',
                'center': 'ttcatcatcattcaacggtg',
                'window': 92,
            },
            'exon4_150PE_2': {
                'refseq': 'catcaccacgctccattatccagccccaaagcgcaacaagcccactgtctatggtgtgtcccccaactacgacaagtgggagatggaacgcacggacatcaccatgaagcacaagctgggcgggggccagtacggggaggtgtacgagggcgtgtggaagaaatacagcctgacggtggccgtgaagaccttgaaggtaggctgggactgccgggggtgcccagggtacgtggggcaag',
                'center': 'catgaagcacaagctgggcg',
                'window': 119,
            },
            'exon5': {
                'refseq': 'gcgctgaagctccattttgcattaactagtcaagtacttacccactgaaaagcacttcctgaaataatttcaccttcgtttttttccttctgcaggaggacaccatggaggtggaagagttcttgaaagaagctgcagtcatgaaagagatcaaacaccctaacctggtgcagctccttggtgagtaagcccggggctctgaagagagggtctcgc',
                'center': 'cttgaaagaagctgcagtca',
                'window': 55,
            },
            'exon6': {
                'refseq': 'ccacgtgttgaagtcctcgttgtcttgttggcaggggtctgcacccgggagcccccgttctatatcatcactgagttcatgacctacgggaacctcctggactacctgagggagtgcaaccggcaggaggtgaacgccgtggtgctgctgtacatggccactcagatctcgtcagccatggagtacctggagaagaaaaacttcatccacaggtaggggcctggccaggcagcctgcgccatggagtcacagggcgtgg',
                'center': 'accggcaggaggtgaacgcc',
                'window': 110,
            },
            'exon7': {
                'refseq': 'ggaaggttggccaggagctctcatgggtgaacattttcctttcttagagatcttgctgcccgaaactgcctggtaggggagaaccacttggtgaaggtagctgattttggcctgagcaggttgatgacaggggacacctacacagcccatgctggagccaagttccccatcaaatggactgcacccgagagcctggcctacaacaagttctccatcaagtccgacgtctggggtaagggctgctgctgcactgaagtggtccttcctg',
                'center': 'caggttgatgacaggggaca',
                'window': 114,
            },
            'exon8': {
                'refseq': 'gtgaaatgctacacatcttgaacagcctttctctttcggttttctttcagcatttggagtattgctttgggaaattgctacctatggcatgtccccttacccgggaattgacctgtcccaggtgtatgagctgctagagaaggactaccgcatggagcgcccagaaggctgcccagagaaggtctatgaactcatgcgagcatgtaagccttcctcagcctgttctcacgagtatatgtgggcattcc',
                'center': 'gacctgtcccaggtgtatga',
                'window': 98,
            },
            'exon9': {
                'refseq': 'agccccgtattgctagccagatctcatggatgatctgacttgggtttcatctgtccaggttggcagtggaatccctctgaccggccctcctttgctgaaatccaccaagcctttgaaacaatgttccaggaatccagtatctcagacggtaaagtacccatcccggggtacctgca',
                'center': 'ccctctgaccggccctcctt',
                'window': 68,
            },
            'invivo_exon4': {
                'refseq': 'catcaccacgctccattatccagccccaaagcgcaacaagcccactgtctatggtgtgtcccccaactacgacaagtgggagatggaacgcacggacatcaccatgaagcacaagctgggcgggggccagtacggggaggtgtacgagggcgtgtggaagaaatacagcctgacggtggccgtgaagaccttgaaggtaggctgggactgccgggggtgcccagggtacgtggggcaag',
                'center': 'catgaagcacaagctgggcg',
                'window': 119,
            },
            'invivo_exon9': {
                'refseq': 'agccccgtattgctagccagatctcatggatgatctgacttgggtttcatctgtccaggttggcagtggaatccctctgaccggccctcctttgctgaaatccaccaagcctttgaaacaatgttccaggaatccagtatctcagacggtaaagtacccatcccggggtacctgca',
                'center': 'ccctctgaccggccctcctt',
                'window': 68,
            },
        }

        return dict_abl1_info[exon]