# Infection Rate Models for COVID-19: Model Risk and Public Health News Sentiment Exposure Adjustments

This repository contains the code for the paper "Infection Rate Models for COVID-19: Model Risk and Public Health News Sentiment Exposure Adjustments", by Ioannis Chalkiadakis, Kevin Hongxuan Yan, Gareth W. Peters, Pavel V. Shevchenko.

Please cite the paper as follows:

_Chalkiadakis I, Yan H, Peters GW, Shevchenko PV (2021) Infection rate models for COVID-19: Model risk and public health news sentiment exposure adjustments. PLoS ONE 16(6): e0253381. https://doi.org/10.1371/journal.pone.0253381 _

The available sentiment data contain the sentiment entropy index per news source, as well as the global weighted sentiment index that was used in the studies.

The COVID-19 data (infected cases) for the studies are also available and taken from https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series.

To obtain the sentiment index used in the studies, run the script 'reporting_weighted_sentiment_index.py' in covid19modelrisk/python/.

To obtain the MCMC results for each model run the corresponding .r file in covid19modelrisk/R/.

Please check the file paths in the code and modify according to your system if necessary.

Tested with Python v3.6.7 and R v4.0.








Copyright @ Authors March 2021

