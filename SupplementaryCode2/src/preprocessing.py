import subprocess, os, datetime, gzip
import pandas as pd
import numpy as np

from tqdm import tqdm
from Bio import SeqIO

class Preprocess:
    def __init__(self, data_path:str, data_format:str='fq.gz'):
        '''Processing NGS raw data for analysis'''

        self.raw_data  = data_path
        self.processed = data_path
        self.data_fmt  = data_format

        self.file_name = data_path.split('/')[-1].replace(f'.{data_format}', '')
        
        self.temp_file = []


    def _check_requirements(self,) -> None:
        '''Check requirements and dependencies for pipeline'''

        if self.data_fmt not in ['fq.gz', 'fastq.gz', 'fq', 'fastq']:
            raise ValueError('Not available format. Please check your input file.')



    def trim(self, finder:str, error:int=0) -> str:
        '''step1 fastq.gz 파일을 cutadapt로 'AAAAAATTCTAG' 찾아서 trimming (3') > .fastq 파일로 저장
        Trimming이 된 후 저장된 파일의 경로 이름을 return 한다. '''

        trimmed = self.processed.replace(f'.{self.data_fmt}', f'_trimmed.{self.data_fmt}')

        command = f'cutadapt -a {finder} -o {trimmed} {self.processed} -e {error}'
        subprocess.call(command, shell=True)

        self.processed = trimmed
        self.temp_file.append(trimmed)

        return trimmed
    
    
    def revcom(self, ) -> str:
        '''step2 seqkit을 이용해서 revcom read로 만들어주기'''


        revcom_file = self.processed.replace(f'.{self.data_fmt}', f'_revcom.{self.data_fmt}')

        command = f'seqkit seq --seq-type DNA -r -p {self.processed} -o {revcom_file}'
        subprocess.call(command, shell=True)

        self.processed = revcom_file
        self.temp_file.append(revcom_file)

        return revcom_file
    

    def to_fasta(self, gzip=True) -> str:
        '''convert fastq to fasta'''

        if   gzip == True : fa_fmt = 'fa.gz'
        elif gzip == False: fa_fmt = 'fa'

        fq = self.processed
        fa = self.processed.replace(self.data_fmt, fa_fmt)

        command = f'seqkit fq2fa {fq} -o {fa}'
        subprocess.call(command, shell=True)

        self.processed = fa
        self.temp_file.append(fa)

        self.data_fmt = fa_fmt

        return fa
    

    
    def finalize(self, save_path:str=None) -> str:
        '''Clean-up temp files and rename final processed data file.'''

        # Delete temp files
        for f in self.temp_file[:-1]:
            if os.path.isfile(f): os.remove(f)

        # Rename final file and move
        final_file_name = f'pp_{self.file_name}.{self.data_fmt}'
        
        if save_path == None:
            command = f'mv {self.processed} {final_file_name}'
        
        else:
            command = f'mv {self.processed} {save_path}/{final_file_name}'

        subprocess.call(command, shell=True)
        
        file_path = command.split(' ')[-1]

        return file_path
    


def make_df_umi(list_barcode:list, data_path:str, len_umi:int) -> pd.DataFrame:
    """NGS read file에서 barcode별로 umi를 구분하고, 읽힌 수를 정리한
    DataFrame을 만들어주는 함수

    Args:
        list_barcode (list): List containing barcodes. pd.Series also acceptable.
        data_path (str): The path of NGS data file. FASTQ or FASTA file can be used.
        len_umi (int): The length of UMI for counting.

    Raises:
        ValueError: NGS data format or path error. 
        ValueError: The lengths of barcode error. Barcode length should be identical.
        ValueError: No barcode error. Check your barcode list.

    Returns:
        _type_: pd.DataFrame
    """    
     
    # Input checker: Check and Determine data file format
    def _check_input(data_path:str):
        '''Input checker: Check and Determine data file format'''

        if   data_path.endswith('.fq.gz')   : data_format = 'fastq'; isgzip = True
        elif data_path.endswith('.fa.gz')   : data_format = 'fasta'; isgzip = True
        elif data_path.endswith('.fastq.gz'): data_format = 'fastq'; isgzip = True
        elif data_path.endswith('.fasta.gz'): data_format = 'fasta'; isgzip = True
        
        elif data_path.endswith('.fq')      : data_format = 'fastq'; isgzip = False
        elif data_path.endswith('.fa')      : data_format = 'fasta'; isgzip = False
        elif data_path.endswith('.fastq')   : data_format = 'fastq'; isgzip = False
        elif data_path.endswith('.fasta')   : data_format = 'fasta'; isgzip = False

        else: raise ValueError('Not supported data format. Please check your input: data_path')

        return data_format, isgzip
    
    data_format, isgzip = _check_input(data_path)

    # Input checker: Check barcode length. The length of barcodes should be identical.
    list_bc_len = [len(bc) for bc in list_barcode]
    if np.std(list_bc_len) != 0: raise ValueError('Please check your input: The lengths of barcode is not identical')
    if len(list_barcode)   == 0: raise ValueError('Please check your input: No barcde found in list_barcode')
    len_bc = list_bc_len[0]
    

    # Step1: Make dictionary containing Barcodes and founded UMIs
    dict_bc = {}
    for bc in list_barcode: dict_bc[bc] = {}

    if isgzip: 
        with gzip.open(data_path, 'rt') as handle:
            list_seq = [str(s.seq) for s in SeqIO.parse(handle, data_format)]

    else:
        list_seq = [str(s.seq) for s in SeqIO.parse(data_path, data_format)]

    for _seq in tqdm(list_seq,
                total = len(list_seq),        ## 전체 진행수
                desc = 'Barcode/UMI sorting', ## 진행률 앞쪽 출력 문장
                ncols = 70,                   ## 진행률 출력 폭 조절
                ascii = ' =',                 ## 바 모양, 첫 번째 문자는 공백이어야 작동
                leave = True
                ):
        
        _bc  = _seq[:len_bc]
        _umi = _seq[-len_umi:]
        
        if _bc in dict_bc: 
            if _umi in dict_bc[_bc]: dict_bc[_bc][_umi] += 1
            else                   : dict_bc[_bc][_umi] = 1
        
        else: continue

    # Step2: Make DataFrame as output
    list_df_temp = []

    for bc in tqdm(dict_bc,
                total = len(dict_bc),  ## 전체 진행수
                desc = 'Make output ', ## 진행률 앞쪽 출력 문장
                ncols = 70,            ## 진행률 출력 폭 조절
                ascii = ' =',          ## 바 모양, 첫 번째 문자는 공백이어야 작동
                leave = True
                ):
        
        list_bc  = []
        list_umi = []
        list_cnt = []
        
        for umi in dict_bc[bc]:
            list_bc.append(bc)
            list_umi.append(umi)
            list_cnt.append(dict_bc[bc][umi])
            
        df_temp = pd.DataFrame()
        df_temp['Barcode'] = list_bc
        df_temp['UMI']     = list_umi
        df_temp['count']   = list_cnt
        
        list_df_temp.append(df_temp)
        
        
    df_out = pd.concat(list_df_temp).reset_index(drop=True)

    return df_out



class MAGeCKanalyzer:
    def __init__(self, ):
        '''Check dependencies and initializing'''

        # Requirements
        # MAGeCK, data directory structures
        # logging?? 


    def setup(self, lib_reference:str, dmso_umi_path:str, tki_umi_path:str, var_type:str='AA_var', rpm:bool=True):
        '''Setup for MAGeCK analysis'''

        self.df_dmso  = self._make_template(lib_reference, dmso_umi_path, var_type)
        self.df_count = self._make_mageck_input(self.df_dmso, data=tki_umi_path, rpm=rpm)

        return self.df_count

    def mageck(self, input_file:str, name:str, control:str='control', test:str='test', save_path:str=None):

        if save_path == None:
            save_dir = f'mageck_result_{name}'
        else:
            save_dir = f'{save_path}/mageck_result_{name}'
        
        os.makedirs(save_dir, exist_ok=True)

        command = f'mageck test -k {input_file} -t {test} -c {control} -n {save_dir}/{name}'

        subprocess.call(command, shell=True)

        df_mageck_result = self._mageck2df(name=name, save_dir=save_dir)
        df_mageck_result.to_csv(f'{save_dir}/{name}_summary.csv')

        return df_mageck_result
    
    def _mageck2df(self, name:str, save_dir:str) -> pd.DataFrame:
        '''MAGeCK results 파일을 DataFrame으로 만들어서 return'''

        
        result = pd.read_csv(f'{save_dir}/{name}.gene_summary.txt', sep='\t').set_index('id')

        pos  = result['pos|score']
        neg  = result['neg|score']
        
        result['p-value'] = [min(pos.loc[i], neg.loc[i]) for i in result.index]

        return result


    def _make_template(self, lib_reference:str, dmso_umi_path:str, var_type:str='AA_var') -> pd.DataFrame:

        if var_type not in ['SNV_var', 'AA_var']:
            raise ValueError('Not available var_type. Please select SNV_var or AA_var')

        df_ref = pd.read_csv(lib_reference).set_index('Barcode')
        df_id = df_ref[[var_type]]
        df_id.columns = ['Gene']

        # DMSO control dataframe
        df_dmso = pd.read_csv(dmso_umi_path)
        df_dmso.columns = ['Barcode', 'UMI_dedup', 'control']

        # SNV sum: DMSO control
        list_gene = [df_id['Gene'].loc[bc] for bc in df_dmso['Barcode']]
        df_dmso.insert(0, 'Gene', list_gene)

        df_dmso['Barcode-UMI'] = df_dmso['Barcode'] + '_' + df_dmso['UMI_dedup']
        df_dmso = df_dmso.set_index('Barcode-UMI')
        df_dmso = df_dmso.drop(['Barcode', 'UMI_dedup'], axis=1)

        return df_dmso


    def _make_mageck_input(self, df_dmso:pd.DataFrame, data:str, rpm:bool=True) -> pd.DataFrame: 
        '''Read count and make input file for MAGeCK'''

        # read sample file
        df_umi = pd.read_csv(data)
        df_umi['Barcode-UMI'] = df_umi['Barcode'] + '_' + df_umi['UMI_dedup']
        
        df_umi = df_umi.set_index('Barcode-UMI')
        df_umi = df_umi.drop(['Barcode', 'UMI_dedup'], axis=1)
        df_umi.columns = ['test']

        # make SNV DataFrame
        df_raw = pd.concat([df_dmso, df_umi], axis=1, join='inner')

        if rpm == True:
            # make RPM file
            df_rpm = df_raw.copy()
            df_rpm['test']    = df_rpm['test'] * 1000000 / np.sum(df_rpm['test'])
            df_rpm['control'] = df_rpm['control'] * 1000000 / np.sum(df_rpm['control'])
            
            return df_rpm
        
        else:
            return df_raw




def pp_log(func):
    def wrapper(*args, **kwargs):
        pass
        


