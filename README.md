# CML VUS screening

본 코드는 CML VUS screening 실험 결과를 분석하기 위한 code를 정리한 것이다. 

### SuppleCode1: epegRNA abundance library design  
pegRNA abundance-based CML VUS screening을 위한 library를 디자인하는 코드이다. Python package [GenET](https://github.com/Goosang-Yu/genet)에 구현된 [DeepPrime-FT (Yu et al., Cell, 2023)](https://doi.org/10.1016/j.cell.2023.03.034)를 활용하여 각 variants를 만들 수 있는 가능한 모든 pegRNA를 디자인하고 그들의 prime editing efficiency를 예측한다. 그리고 최적의 pegRNA들을 선정하여 library를 구성하였다. 

### SuppleCode2: epegRNA abundance screening analysis  
pegRNA abundance-based CML VUS screening의 결과로 얻은 데이터를 전처리하고 분석하는데 사용한 코드이다. SRA에 업로드 되어있는 raw FASTQ 파일을 분석할 부분만 남기도록 전처리한 후, barcode sequence마다의 read count를 세고, UMI clustering을 해서 variants 마다의 read count를 정리하였다. 이렇게 얻어낸 variants 별 read count는 MAGeCK pipeline을 사용해 분석하였다. 

### SuppleCode3: SynPrime library design  
SynPrime-based CML VUS screening을 위한 library를 디자인하는 코드이다. 

### SuppleCode4: SypPrime screening analysis  
SynPrime-based CML VUS screening의 결과로 얻은 데이터를 전처리하고 분석하는데 사용한 코드이다. Alignment는 [CRISPResso2](https://github.com/pinellolab/CRISPResso2)를 활용하였으며, alignment 후 만들어진 frequency talbe을 토대로 이후 분석을 진행하였다. 각 variants를 만들기 위한 pegRNA와 일치하는 read와 그 외의 나머지 reads의 패턴도 분석하였다. 또한, screening 결과로 찾아낸 hit variants에 대해서 추가적인 cell level validation에 사용된 분석용 코드도 포함되어 있다. 

### SuppleCode5: In vivo validation analysis  
SynPrime screening으로 밝혀낸 CML TKI resistant variants들을 in vivo xenograft model에서 검증한 실험 결과에 대해 분석한 코드이다. 

### SuppleCode6: Off-target validation  
이번 연구에서 사용한 library 내의 pegRNA 중, off-target effect가 있을 것으로 예상되는 pegRNA들을 선정하고, 이를 분석하는데 사용한 코드. 

# Raw sequencing data
본 분석에 사용된 raw NGS file들은 모두 SRA에 올라가 있다 (BioProject; [PRJNA1048659](https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA1048659)). 각 pipeline에 사용된 FASTQ 파일들의 accession number는 각 SupplementaryCode의 README를 참고하면 된다. 

# Environments
본 코드들은 Ubuntu22.04 LTS 또는 CentOS7의 환경에서 테스트 되었다. 

# Requirements
- Python > 3.6
- biopython
- pandas
- genet
