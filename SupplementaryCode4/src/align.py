import sys, os
import subprocess
import pandas as pd
from glob import glob

class ABL1VUS:
    def __init__(self, sample_id:str, r1, r2, exon:str, out_dir:str):
        '''CML VUS screening에서 ABL1의 variants counting을 위한 함수
        CRISPResso를 이용한 counting을 하고, 그 output을 이용해 variants list별로 정리된 
        output 파일을 생성한다. '''
        
        self._check_input(exon)

        self.info = self._exon_align_info(exon)

        self.refseq = self.info['refseq']
        self.center = self.info['center']
        self.window = self.info['window']

        data    = f'-r1 {r1} -r2 {r2} -n {sample_id}'
        align   = f'-a {self.refseq} -an {exon} -g {self.center} --plot_window_size {self.window}'
        output  = f'-o {out_dir} --file_prefix {sample_id}'
        fixed   = f'--suppress_plots --suppress_report'
        
        command = f'CRISPResso {data} {align} {output} {fixed}'

        subprocess.run([command], shell=True, check=True)



    def _check_input(self, exon):
        if exon not in ['exon4', 'exon5', 'exon6', 'exon7', 'exon8', 'exon9']:
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
        }

        return dict_abl1_info[exon]

base_dir = '/extdata2/CML_vus/10_NBT_revision/1_raw_data/X201SC24040832-Z01-F001/01.RawData/'

list_sample = [

                # 'R1_DMSO_PE2_E5',
                # 'R1_DMSO_PE4_E5',
                # 'R1_DMSO_PE4K_E5',
                # 'R1_DMSO_KCL22_PE4K_E5',
                # 'R1_Ima_PE2_E5',
                # 'R1_Ima_PE4_E5',
                # 'R1_Ima_PE4K_E5',
                # 'R1_Ima_KCL22_PE4K_E5',
                # 'R1_Nilo_PE2_E5',
                # 'R1_Nilo_PE4_E5',
                # 'R1_Nilo_PE4K_E5',
                # 'R1_Nilo_KCL22_PE4K_E5',
                # 'R1_Bosu_PE2_E5',
                # 'R1_Bosu_PE4_E5',
                # 'R1_Bosu_PE4K_E5',
                # 'R1_Bosu_KCL22_PE4K_E5',
                # 'R1_Dasa_PE2_E5',
                # 'R1_Dasa_PE4_E5',
                # 'R1_Dasa_PE4K_E5',
                # 'R1_Dasa_KCL22_PE4K_E5',
                # 'R1_Pona_PE2_E5',
                # 'R1_Pona_PE4_E5',
                # 'R1_Pona_PE4K_E5',
                # 'R1_Pona_KCL22_PE4K_E5',
                # 'R1_Asci_PE2_E5',
                # 'R1_Asci_PE4_E5',
                # 'R1_Asci_PE4K_E5',
                # 'R1_Asci_KCL22_PE4K_E5',
                # 'R1_Ima0_5x_PE4K_E5',
                # 'R1_Ima2x_PE4K_E5',
                # 'R1_Asci0_5x_PE4K_E5',
                # 'R1_Asci2x_PE4K_E5',
                # 'PE4K_Unedit_E5',

                'R1_DMSO_PE4K_E5_17cycle',

                # 'R1_DMSO_PE4K_E5_20cycle',
                # 'R1_DMSO_PE4K_E5_23cycle',


            ]


for sample_id in list_sample:
    
    fq_path = f'{base_dir}/{sample_id}'
    files = list(glob(f'{fq_path}/*.fq.gz'))
    
    exon_num = 5
    
    r1 = files[0]
    r2 = files[1]

    abl_e8 = ABL1VUS(sample_id, r1, r2, exon=f'exon{exon_num}', out_dir='./')