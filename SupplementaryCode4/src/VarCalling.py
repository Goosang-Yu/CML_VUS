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


class ReadPatternAnalyzer:
    def __init__(self,):
        """### Mutation type에 따른 분류
        Substitution / insertion / deletion / complex에 따라서 파일을 나눠서 저장하기

        Case1: Aligned_Sequence에만 '-'가 들어있는 경우 > Deletion  
        Case2: Reference_Sequence에만 '-'가 들어있는 경우 > Insertion  
        Case3: 둘 다 '-'가 들어있지 않은 경우 > Substitution (or WT)  
        Case4: 그 외의 나머지, 여러개의 mutation이 복합적으로 들어있는 것 > Complex
        """        

        pass

    
    def run(self, freq_table:str, ref_info:str) -> pd.DataFrame:
        """CRISPResso2 분석 완료된 파일을 기준으로 read pattern을 분석하는 pipeline 실행.

        Args:
            freq_table (str): CRISPResso2 결과로 나온 frequency table
            ref_info (str): 예상되는 read에 대해서 white list를 만들어둔 것. 

        Returns:
            pd.DataFrame: Read pattern 정보가 추가된 DataFrame
        """        

        df_freq = pd.read_csv(freq_table, sep='\t')
        df_ref = pd.read_csv(ref_info)

        is_ins = df_freq['Reference_Sequence'].str.contains('-')
        is_del = df_freq['Aligned_Sequence'].str.contains('-')

        df_ins = df_freq[is_ins & ~is_del]
        df_del = df_freq[~is_ins & is_del]
        df_sub = df_freq[~is_ins & ~is_del]

        list_idx  = list(df_sub.index) + list(df_ins.index) + list(df_del.index)
        df_complx = df_freq.drop(list_idx)
        df_complx['mut_type'] = 'Complex'
        df_complx['mut_class'] = 'Complex'

        # classification
        df_sub_type = self._classify_substitutions(df_sub, df_ref)
        df_ins_type = self._classify_indel(df_ins, 'insertion')
        df_del_type = self._classify_indel(df_del, 'deletion')

        df_merge = pd.concat([df_sub_type, df_ins_type, df_del_type, df_complx])
        df_merge = df_merge.sort_index()

        df_merge = df_merge[['#Reads', '%Reads', 'mut_type', 'mut_class']]
        
        return df_merge


    def _count_mismatches(self, seq1:str, seq2:str) -> int:
        """Compare 2 sequences and return the number of mismatches"""

        cnt = 0

        for a, b in zip(list(seq1), list(seq2)):
            if a != b: cnt += 1
        
        return cnt


    def _classify_substitutions(self, df_reads:pd.DataFrame, df_ref:pd.DataFrame):
        """Substitution만 있는 것에서, WT과 SynPrime이 정확히 일어난 것을 분리하고,
        그 외에 unintended edit만 포함된 것을 따로 DataFrame으로 분류하기.

        Args:
            df_reads (pd.DataFrame): _description_
            df_ref (pd.DataFrame): _description_

        Returns:
            _type_: _description_
        """        

        # Step1: make dictionary for read count
        dict_ref = {}

        for idx in df_ref.index:
            dict_ref[df_ref.iloc[idx].RefSeq] = df_ref.iloc[idx].Label

        # Step2: classify substitution types
        list_sub_type = []
        list_sub_class= []

        for idx in tqdm(df_reads.index, total=len(df_reads.index), 
                        desc='Classify pattern', ncols=100, leave=False):

            data = df_reads.loc[idx]

            try: 
                mut_type = dict_ref[data.Aligned_Sequence]
                list_sub_type.append(mut_type)

                if mut_type in ['Intended_only', 'Synony_only']: 
                    list_sub_class.append('Single_edit')
                else:
                    list_sub_class.append(dict_ref[data.Aligned_Sequence])
            
            # SynPrime으로 생길 수 없는 product에 대해서 분류
            except:
                cnt = self._count_mismatches(data.Aligned_Sequence, data.Reference_Sequence)
                list_sub_type.append(f'sub{cnt}')
                
                if cnt > 4: list_sub_class.append(f'sub5more')
                else      : list_sub_class.append(f'sub{cnt}')


        df_reads_type = df_reads.copy()
        df_reads_type['mut_type'] = list_sub_type
        df_reads_type['mut_class'] = list_sub_class

        return df_reads_type
        
    def _find_indel_sequence(self, aligned_seq:str, ref_seq:str, type:str='insertion') -> str:
        """Aligned reads의 insertion/deletion pattern을 분석하는 handler.

        Args:
            aligned_seq (str): _description_
            ref_seq (str): _description_
            type (str): Select'insertion' or 'deletion'

        Returns:
            str: Sequence of insertion or deletion.
        """

        if   type == 'insertion': pass
        elif type == 'deletion' :
            aligned_seq, ref_seq = ref_seq, aligned_seq
        else: 
            raise ValueError('Not available type. Select "insertion" or "deletion"')

        list_ins = []
        isInsert = False

        for a, b in zip(list(aligned_seq), list(ref_seq)):

            if b == '-':
                isInsert = True
                list_ins.append(a)

            else:
                if isInsert == False: pass
                else: break

        return ''.join(list_ins)

    def _classify_indel(self, df_reads:pd.DataFrame, type:str):
        """Frequency table read에서 insertion / deletion의 pattern을 분석하는 코드

        Args:
            df_reads (pd.DataFrame): _description_
            type (str): Select'insertion' or 'deletion'

        Returns:
            _type_: _description_
        """        

        # make list for inserted seq
        list_indel_seq = []

        for idx in df_reads.index:

            data = df_reads.loc[idx]
            aligned_seq = data.Aligned_Sequence
            ref_seq     = data.Reference_Sequence

            # Find inserted sequence
            indel_seq = self._find_indel_sequence(aligned_seq, ref_seq, type)
            
            list_indel_seq.append(f'{type[:3]}{len(indel_seq)}:'+indel_seq)

        df_out = df_reads.copy()
        df_out['mut_type'] = list_indel_seq
        df_out['mut_class'] = type[:3]

        return df_out


def single_clone_var_freq(sample_id:str, freq_table:str, wt_seq:str, edit_seq:str, intended_only:str) -> pd.DataFrame:
    
    wt_seq   = wt_seq.upper()
    edit_seq = edit_seq.upper()
    
    # Step1: read CRISPResso aligned & reference file
    df = pd.read_csv(freq_table, sep = '\t')
    
    sample_name = os.path.basename(freq_table).replace('.txt', '')
    dict_out = {wt_seq: 0, edit_seq: 0, intended_only:0, 'Others': 0}

    print(f'[Info] Start - {sample_name}')

    # Step2: read count
    for i in df.index:
        data = df.iloc[i]
        seq  = data['Aligned_Sequence']
        cnt  = int(data['#Reads'])
        
        try   : dict_out[seq] += cnt
        except: dict_out['Others'] += cnt

    dict_out['WT'] = dict_out.pop(wt_seq)
    dict_out['Variant'] = dict_out.pop(edit_seq)
    
    try   : dict_out['Intended_only'] = dict_out.pop(intended_only)
    except: dict_out['Intended_only'] = 0

    df_out = pd.DataFrame.from_dict(dict_out, orient='index').T
    df_out = df_out.rename(index={0: sample_id})

    return df_out