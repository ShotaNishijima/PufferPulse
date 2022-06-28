# Code Share of Nishijima et al. (202X)  

This repository shares the following programming codes used in Nishijima et al. (202X)
 
- surfnet_peak_model.cpp: TMB code for modeling pulse dynamics of catch rate with the Guassian function
- surnet_peak_modeling.R: Source code for various functions
- 01_pulse_modeling.r: Running code to model the pulse dynamics of catch rate for standardized abundance index
- 02_vpa_hindcasting.r: Running code to analyze VPA with short-term forecasting and hindcasting

## Data share

This repository also shares the datasets used in Nishijima et al. (202X)

- data/surnet_survey2021.csv: Data of surfnet survey from 2004 to 2021
- data/caa.csv: Catch-at-age data of the Ise-Mikawa Bay stock of Japanese pufferfish (<i>Takifugu rubripes</i>)
- data/waa.csv: Weight-at-age data of this stock
- data/maa.csv: Maturity-at-age data of this stock
- data/M.csv: Natural mortality-at-age data of this stock
- data/cpue_age1: CPUE data for age-1 fish of this stock
- data/Recruit_estimate2020: Estimates of hathery-reard reruitment of this stock

## Results

The folder "res" includes all results of figures and tables obtained by running "01_pulse_modeling.r" and "02_vpa_hindcasting.r".

## Citation 

### Code and surf-net survey data
- Nishijima S., Suzuki S., Fukuta R., and Okada M. (2022) Modeling pulse dynamics of juvenile fish enables the short-term forecasting of population dynamics in Japanese pufferfish: A latent variable approach. bioRxiv [doi: https://doi.org/10.1101/2022.01.26.477932]

### Data for stock assessment (caa, waa, maa, M, cpue_age1, Recruit_estimate2020)
- Suzuki, S., Yamamoto, T., Sawayama, S., and Nishijima, S. (2021) Stock assessment and evaluation for the Ise-Mikawa Bay population of ocellate puffer (fiscal year 2020). In Marine Fisheries Stock Assessment and Evaluation for Japanese Waters. Tokyo, Japan (in Japanese). pp. 1–44. Available from http://abchan.fra.go.jp/digests2020/index.html

### VPA analysis with the "frasyr" package
- Please run 'citation("frasyr")' and see https://github.com/ichimomo/frasyr

### About age-1 cpue index 
- Nishijima, S., Suzuki, S., Ichinokawa, M., and Okamura, H. (2019) Integrated multi-timescale modeling untangles anthropogenic, environmental, and biological effects on catchability. Canadian Journal of Fisheries and Aquatic Sciences 76: 2045–2056. [doi: https://doi.org/10.1139/cjfas-2018-0114]
