# Infection Rate Models for COVID-19: Model Risk and Public Health News Sentiment Exposure Adjustments

This repository contains the code for the paper "Infection Rate Models for COVID-19: Model Risk and Public Health News Sentiment Exposure Adjustments", by Ioannis Chalkiadakis, Kevin Hongxuan Yan, Gareth W. Peters, Pavel V. Shevchenko.


The available sentiment data contain the sentiment entropy index per news source, as well as the global weighted sentiment index. 

The infected COVID-19 data for the studies are also available and taken from https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series.

To obtain the sentiment index used in the studies, run the script 'reporting_weighted_sentiment_index.py' in covid19modelrisk/python/.

To obtain the MCMC results for each model run the corresponding .r file in covid19modelrisk/R/.

Please check the file paths in the code and modify according to your system if necessary.

Tested with Python v3.6.7 and R v4.0.








Copyright @ Authors March 2021

