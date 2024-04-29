import sys, os
import pandas as pd

from tqdm import tqdm
from glob import glob
from scipy.stats import fisher_exact

def make_count_file(freq_table:str, var_ref:str) -> pd.DataFrame:
    """CRISPResso2를 이용해서 각 variants마다의 read를 alignment 한 파일에서 read count를 가져온다. 
    
    Args:
        freq_table (str): CRISPResso로 생성된 frequency table 파일의 경로
        var_ref (str): 분석하고자 하는 variant가 포함된 서열과 해당 variants의 정보가 담긴 reference file의 경로

    Returns:
        pd.DataFrame: _description_
    """    
    
    # Step1: read CRISPResso aligned & reference file
    df_ref = pd.read_csv(var_ref)
    df = pd.read_csv(freq_table, sep = '\t')
    
    sample_name = os.path.basename(freq_table).replace('.txt', '')
    dict_out = {}

    for ref_seq in df_ref['RefSeq']:
        dict_out[ref_seq]  = 0

    print(f'\n[Info] Start - {sample_name}')
    print('[Info] Length of variants:', len(dict_out))

    # Step2: read count
    for i in tqdm(df.index, desc=f'[Info] Read counting: {sample_name}'):
        data = df.iloc[i]
        seq  = data['Aligned_Sequence']
        cnt  = int(data['#Reads'])
        
        try   : dict_out[seq] += cnt
        except: dict_out['No_matched'] += cnt

    # Step3: make output
    list_count = []

    for ref_seq in df_ref['RefSeq']:
        list_count.append(dict_out[ref_seq])

    df_out = df_ref.copy()
    df_out['count'] = list_count

    return df_out


def read_statistics(var_sample:str, background:str) -> pd.DataFrame:
    """_summary_

    Args:
        var_sample (str): _description_
        background (str): _description_

    Raises:
        ValueError: _description_

    Returns:
        pd.DataFrame: _description_
    """    

    ## Check var_sample과 backgound의 variants list가 완전히 동일한지 확인!
    if False in list(df_UE['RefSeq'] == df_sample['RefSeq']):
        raise ValueError('Not matched between sample and background. Please check your input files.')

    # Step1: Unedit dictionary 만들기

    df_UE = pd.read_csv(background)
    
    UE_WT_read  = df_UE[df_UE['Label']=='WT_refseq']['count'].iloc[0]
    UE_SynPE_df = df_UE[df_UE['Label']=='SynPE'].reset_index(drop=True)

    UE_SynPE_dict = {}

    for i in UE_SynPE_df.index:
        UE_data = UE_SynPE_df.iloc[i]
        UE_SynPE_dict[UE_data['RefSeq']] = int(UE_data['count'])

    # Step2: 각 Stat file마다 odds / p-value column 추가하기

    f_name = os.path.basename(var_sample).replace('.csv', '')
    print('Analysis:', f_name)
    
    df_sample = pd.read_csv(var_sample)
    df_synpe  = df_sample[df_sample['Label']=='SynPE'].reset_index(drop=True).copy()

    total_cnt_wtseq = df_sample[df_sample['Label']=='WT_refseq']['count'].iloc[0]
    total_cnt_edseq = df_synpe['count'].sum()

    list_sample_WT  = [] # list_sample_wt_cnt
    list_sample_RPM = []
    list_UE_SynPE   = []
    list_UE_WT_read = []
    list_odds_ratio = []
    list_pvalue     = []

    for i in df_synpe.index:

        data = df_synpe.loc[i]

        sample_SynPE_cnt = int(df_synpe.loc[i]['count'])
        sample_SynPE_rpm = sample_SynPE_cnt*1000000/total_cnt_edseq
        unedit_SynPE_cnt = UE_SynPE_dict[data['RefSeq']]

        # Odds ratio 계산
        odds = ((sample_SynPE_cnt+1)/(total_cnt_wtseq+1))/((unedit_SynPE_cnt+1)/(UE_WT_read+1))
        
        # Fisher test p-value 계산
        fisher_df = pd.DataFrame({'intended_edit':[sample_SynPE_cnt, unedit_SynPE_cnt], 'WT':[total_cnt_wtseq, UE_WT_read]})
        oddsr, p = fisher_exact(table = fisher_df.to_numpy(), alternative = 'two-sided')

        # 각각 계산한 값들을 list에 담기
        list_sample_WT.append(total_cnt_wtseq)
        list_sample_RPM.append(sample_SynPE_rpm)
        list_UE_SynPE.append(unedit_SynPE_cnt)
        list_UE_WT_read.append(UE_WT_read)
        list_odds_ratio.append(odds)
        list_pvalue.append(p)

    df_synpe['Edited_WT_count'] = list_sample_WT
    df_synpe['RPM']             = list_sample_RPM
    df_synpe['UE_SynPE_count']  = list_UE_SynPE
    df_synpe['UE_WT_count']     = list_UE_WT_read
    df_synpe['OR']              = list_odds_ratio
    df_synpe['pvalue']          = list_pvalue

    return df_synpe



