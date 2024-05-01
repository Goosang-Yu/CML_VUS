import joypy
import pandas as pd
import matplotlib.pyplot as plt

def get_class_list(list_var:list, class_file:str, tki:str) -> list:
    """In vivo 데이터 분석을 할 때, 각 variants들의 classification을 가져오기 위한 함수

    Args:
        list_var (list): Class를 알고 싶은 variants들의 list
        class_file (str): 각 variants마다의 classficiation이 정리된 파일 경로. variants_info에 들어있다.
        tki (str): 보려는 TKI의 종류. CML에 대한 6개의 FDA 승인된 약물로 한정.

    Raises:
        ValueError: TKI가 잘못된 것으로 지정되면 error 발생.

    Returns:
        list: 입력해준 variants들에 대한 class들이 담긴 list
    """    

    if tki not in ['Imatinib', 'Nilotinib', 'Bosutinib', 'Dasatinib', 'Ponatinib', 'Asciminib']:
        raise ValueError('Not available TKI. Please check your input.')
    
    class_info = pd.read_csv(class_file).set_index('variant')[tki]
    
    return [class_info[var] for var in list_var]

    

def make_invivo_data(sample:str, tki:str, data_dir:str='data/5_LFC', class_file:str='variants_info/variants_class.csv') -> pd.DataFrame:
    """In vivo data를 classification까지 붙여서 정리한 DataFrame으로 만들어주는 함수.
    여기서 만들어진 DataFrame을 이용해서 joyplot을 그린다. 

    Args:
        sample (str): In vivo sample ID
        tki (str): 사용한 TKI 종류. Classification 정보를 불러올 때 사용된다. 
        data_dir (str, optional): _description_. Defaults to 'data/5_LFC'.
        class_file (str, optional): _description_. Defaults to 'variants_info/variants_class.csv'.

    Returns:
        pd.DataFrame: In vivo data를 classification까지 붙여서 정리한 DataFrame
    """

    data = pd.read_csv(f'{data_dir}/{sample}.csv').set_index('SAAV')
    data['Class'] = get_class_list(list(data.index), class_file, tki)

    data = data.replace({'Resistant':'1_Resistant', 'Intermediate':'2_Intermediate', 'Sensitive': '3_Sensitive'})

    return data[['SNV', 'pos', 'LFC', 'Class']]


def make_invivo_joyplot(data:pd.DataFrame, save_path:str=None, x_range=[-5, 10]):
    """In vivo data를 이용해서 LFC distribution을 볼 수 있는 joyplot을 만드는 함수

    Args:
        data (pd.DataFrame): LFC와 variants 등의 정보가 담겨있는 DataFrame
        save_path (str, optional): 그려진 joyplot을 저장할 경로. Defaults to None.
    """

    plt.figure(dpi= 1200)

    fig, axes = joypy.joyplot(data, column=['LFC'], by="Class", ylim='own', x_range=x_range, figsize=(4,2.3), 
                            color=['#003f5c', '#ef5675', '#ffa600'], legend=False, alpha=0.7, overlap=0.4)

    # Decoration
    plt.rc("font", size=10)
    plt.xlabel('Log Fold Change',  fontsize=10, color='black', alpha=1)
    plt.ylabel('Density', fontsize=10,  color='black', alpha=0.8)

    # plt.show
    if save_path != None:
        plt.savefig(save_path, dpi= 1200)
        plt.close()


