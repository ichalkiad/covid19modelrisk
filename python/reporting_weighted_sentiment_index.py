###############################################################################
# "Infection Rate Models for COVID-19: 
#  Model Risk and Public Health News Sentiment Exposure Adjustments"
#   
#  Ioannis Chalkiadakis, Kevin Hongxuan Yan, Gareth W. Peters, Pavel V. Shevchenko
#
#  Ioannis Chalkiadakis, ic14@hw.ac.uk
#  March 2021
###############################################################################

import json
import time
import logging
import datetime
from .text_utils import get_timeseries_summary_per_day, get_timestamp_gaps_dates
import pickle
import os
from scipy.io import savemat
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy import ndimage
import pandas as pd


if __name__ == "__main__":

	DIR = "./covid19modelrisk/data/sentiment"
	DIR_out_base = "/covid19modelrisk/output/"

	interpolation = True
	# 1-week window for median filter
	median_filter_size = 7
	daily_summary = "IQR"
	kwargs_list = ["nyt", "ecdc", "uscdc", "who", "unece"]
	sentiment_type = ["total"]
	ts = ["ts_token_entropy.pickle"]

	DIR_out_base = DIR_out_base + daily_summary + "/"
	if not os.path.exists(DIR_out_base):
	    os.makedirs(DIR_out_base)

	matlab_output = dict()
	matlab_output_fwd = dict()

	logging.basicConfig(filename=DIR_out_base + "/run_join_ts.log", level=logging.INFO)
	t0 = time.time()

	for senttype in sentiment_type:
	    DIR_out = DIR_out_base + "/{}/".format(senttype)
	    if not os.path.exists(DIR_out):
		os.makedirs(DIR_out)

	    for t in ts:
		data_meta = []
		data_series = []
		# unique date keys
		idxs = []
		for top in range(len(kwargs_list)):
		    series = []
		    extra = []
		    elem = []

		    source = kwargs_list[top]
		    ts_name = "{}/{}_{}/".format(DIR, source, senttype
		    try:
		        with open(ts_name + "/timeseries_elements.pickle", "rb") as f:
		            timeseries_elements = pickle.load(f)
		            metadata = pickle.load(f)
		    except FileNotFoundError:
		        logging.info("{} was not found - check ts construction script".format(ts_name
		                                                                              + "/timeseries_elements.pickle"))
		    meta = []
		    ts_elements = []
		    try:
		        with open(ts_name + "/" + t, "rb") as f:
		            ts_data = pickle.load(f)
		            if len(ts_elements) > 0:
		                telements = ts_elements
		                metadata = meta
		            else:
		                telements = timeseries_elements
		    except FileNotFoundError:
		        logging.info("{} was not found.".format(ts_name+"/"+t))
		        continue
		    series.extend(ts_data)
		    extra.extend(metadata)
		    elem.extend(telements)

		    summary_ts, summary_meta, time_batches = get_timeseries_summary_per_day(series, elem, extra,
		                                                                            summary=daily_summary)
		    idxs.extend([sk for sk in summary_ts.keys() if sk not in idxs])
		    data_series.append(summary_ts)
		    data_meta.append(summary_meta)

		contents = []
		dtype = []
		headers_columns = ["Date", "25%", "50%", "75%"]
		for h in headers_columns:
		    dtype.append((h, (np.str_, 1000)))

		# sorting and weighting
		total_summaries = []
		total_meta = []
		nyt_weights = []
		ecdc_weights = []
		uscdc_weights = []
		who_weights = []
		unece_weights = []
		idxs = sorted(idxs)
		for j in range(len(idxs)):
		    i = idxs[j]
		    # add date
		    row = [str(i)]
		    lower_p = ""
		    median_p = ""
		    upper_p = ""
		    total_sum = 0
		    try:
		        # NYT val and ngram number
		        nyt_val = data_series[0][i]
		        nyt_num = data_meta[0][i]["number of ngrams"]
		        nyt_sum = data_meta[0][i]
		        lower_p += data_meta[0][i]["25%"] + ","
		        median_p += data_meta[0][i]["50%"] + ","
		        upper_p += data_meta[0][i]["75%"] + ","
		    except:
		        nyt_val = 0
		        nyt_num = 0
		        nyt_sum = []
		    try:
		        # ECDC
		        ecdc_val = data_series[1][i]
		        ecdc_num = data_meta[1][i]["number of ngrams"]
		        ecdc_sum = data_meta[1][i]
		        lower_p += data_meta[1][i]["25%"] + ","
		        median_p += data_meta[1][i]["50%"] + ","
		        upper_p += data_meta[1][i]["75%"] + ","
		    except:
		        ecdc_val = 0
		        ecdc_num = 0
		        ecdc_sum = []
		    try:
		        # USCDC
		        uscdc_val = data_series[2][i]
		        uscdc_num = data_meta[2][i]["number of ngrams"]
		        uscdc_sum = data_meta[2][i]
		        lower_p += data_meta[2][i]["25%"] + ","
		        median_p += data_meta[2][i]["50%"] + ","
		        upper_p += data_meta[2][i]["75%"] + ","
		    except:
		        uscdc_val = 0
		        uscdc_num = 0
		        uscdc_sum = []
		    try:
		        # WHO
		        who_val = data_series[3][i]
		        who_num = data_meta[3][i]["number of ngrams"]
		        who_sum = data_meta[3][i]
		        lower_p += data_meta[3][i]["25%"] + ","
		        median_p += data_meta[3][i]["50%"] + ","
		        upper_p += data_meta[3][i]["75%"] + ","
		    except:
		        who_val = 0
		        who_num = 0
		        who_sum = []
		    try:
		        # UNECE
		        unece_val = data_series[4][i]
		        unece_num = data_meta[4][i]["number of ngrams"]
		        unece_sum = data_meta[4][i]
		        lower_p += data_meta[4][i]["25%"] + ","
		        median_p += data_meta[4][i]["50%"] + ","
		        upper_p += data_meta[4][i]["75%"] + ","
		    except:
		        unece_val = 0
		        unece_num = 0
		        unece_sum = []
		    total_sum = nyt_num + ecdc_num + uscdc_num + who_num + unece_num
		    total_summaries.append(nyt_val*(nyt_num/total_sum) + ecdc_val*(ecdc_num/total_sum) +
		                           uscdc_val*(uscdc_num/total_sum) + who_val*(who_num/total_sum) +
		                           unece_val*(unece_num/total_sum))
		    total_meta.append([nyt_sum, ecdc_sum, uscdc_sum, who_sum, unece_sum])
		    nyt_weights.append(nyt_num / total_sum)
		    ecdc_weights.append(ecdc_num / total_sum)
		    uscdc_weights.append(uscdc_num / total_sum)
		    who_weights.append(who_num / total_sum)
		    unece_weights.append(unece_num / total_sum)

		    for ngram_prc in [lower_p, median_p, upper_p]:
		        perc = ngram_prc.replace(",", "-").split('-')[0:-1]
		        per = ""
		        for ngr in range(len(perc)):
		            per += perc[ngr] + "-"
		            # if ngr > 0 and ngr % 10 == 0:
		            #     per += "\n"
		        row.append("{}".format(per[0:-1]))

		    contents.append(tuple(row))

		arr = np.array(contents, dtype=dtype)
		np.savetxt(DIR_out + "ngram_table.csv", arr, delimiter=',', fmt=['%s' for i in range(len(headers_columns))],
		       header=",".join(headers_columns), comments='')

		with open(DIR_out + "/" + t.replace(".pickle", "_joined.pickle"), "wb") as f:
		    pickle.dump(total_summaries, f, pickle.HIGHEST_PROTOCOL)
		    pickle.dump(total_meta, f, pickle.HIGHEST_PROTOCOL)
		    pickle.dump(idxs, f, pickle.HIGHEST_PROTOCOL)

		matlab_output[t.replace(".pickle", "_joined")] = total_summaries
		matlab_output[t.replace(".pickle", "_joined_metadata")] = total_meta
		matlab_output[t.replace(".pickle", "_time")] = [str(i) for i in idxs]

		gap_idx = get_timestamp_gaps_dates(idxs, DIR_out + "/" + t.replace(".pickle", "_gaps.html"))
		gap_indices = [i[0] for i in gap_idx]

		total_summaries_fwd = []
		total_meta_fwd = []
		idxss = []
		nyt_weights_extend = []
		ecdc_weights_extend = []
		uscdc_weights_extend = []
		who_weights_extend = []
		unece_weights_extend = []
		ftxt = open(DIR_out + "/gap_dates.txt", "wt")

		if interpolation:
		    fitting_idx = []
		    interp_idx = []
		    for i in range(len(idxs) - 1):
		        # linear / cubic spline interpolation for missing days
		        diff_days = int((idxs[i + 1] - idxs[i]) / np.timedelta64(1, 'D'))
		        if diff_days == 1:
		            fitting_idx.append(i + len(interp_idx) + 1)
		        else:
		            if len(interp_idx) == 0:
		                interp_idx.extend([i + j for j in range(1, diff_days, 1)])
		            else:
		                interp_idx.extend([fitting_idx[-1] + j for j in range(1, diff_days, 1)])
		            if len(fitting_idx) == 0:
		                fitting_idx.append(i)
		                fitting_idx.append(i + diff_days)
		            else:
		                fitting_idx.append(fitting_idx[-1] + diff_days)

		    f = interpolate.interp1d(fitting_idx, total_summaries, kind="linear")

		for i in range(int((idxs[-1] - idxs[0]) / np.timedelta64(1, 'D')) + 1):
		    if interpolation and i in interp_idx:
		        total_summaries_fwd.append(f(i))
		        total_meta_fwd.append(-1)
		        idxss.append(idxss[-1] + np.timedelta64(1, 'D'))
		        ftxt.write(str(idxss[-1]) + "\n")
		        # for gap days, append zeros when there's a gap
		        nyt_weights_extend.append(0)
		        ecdc_weights_extend.append(0)
		        uscdc_weights_extend.append(0)
		        who_weights_extend.append(0)
		        unece_weights_extend.append(0)
		    elif interpolation and i in fitting_idx:
		        id = fitting_idx.index(i)
		        nyt_weights_extend.append(nyt_weights[id])
		        ecdc_weights_extend.append(ecdc_weights[id])
		        uscdc_weights_extend.append(uscdc_weights[id])
		        who_weights_extend.append(who_weights[id])
		        unece_weights_extend.append(unece_weights[id])
		        total_summaries_fwd.append(total_summaries[id])
		        total_meta_fwd.append(total_meta[id])
		        idxss.append(idxs[id])

		ftxt.close()
		if not os.path.isfile(DIR_out + "/source_weights.pickle"):
		    with open(DIR_out + "/source_weights.pickle", "wb") as f:
		        pickle.dump(nyt_weights_extend, f, pickle.HIGHEST_PROTOCOL)
		        pickle.dump(ecdc_weights_extend, f, pickle.HIGHEST_PROTOCOL)
		        pickle.dump(uscdc_weights_extend, f, pickle.HIGHEST_PROTOCOL)
		        pickle.dump(who_weights_extend, f, pickle.HIGHEST_PROTOCOL)
		        pickle.dump(unece_weights_extend, f, pickle.HIGHEST_PROTOCOL)
		        pickle.dump(idxss, f, pickle.HIGHEST_PROTOCOL)

		# apply median filter to the signal - boundaries are handled by mirroring around the boundary value
		median_filtered_signal_fwd = ndimage.median_filter(total_summaries_fwd, size=((median_filter_size,)),
		                                               mode='mirror')

		with open(DIR_out + "/" + t.replace(".pickle", "_joined_fwd.pickle"), "wb") as f:
		    pickle.dump(total_summaries_fwd, f, pickle.HIGHEST_PROTOCOL)
		    pickle.dump(median_filtered_signal_fwd, f, pickle.HIGHEST_PROTOCOL)
		    pickle.dump(total_meta_fwd, f, pickle.HIGHEST_PROTOCOL)
		    pickle.dump(idxss, f, pickle.HIGHEST_PROTOCOL)

		matlab_output_fwd[t.replace(".pickle", "_joined")] = total_summaries_fwd
		matlab_output_fwd[t.replace(".pickle", "_joined_median_filtered")] = median_filtered_signal_fwd
		matlab_output_fwd[t.replace(".pickle", "_joined_metadata")] = total_meta_fwd
		matlab_output_fwd[t.replace(".pickle", "_time")] = [str(i) for i in idxss]
		matlab_output_fwd[t.replace(".pickle", "_gaps")] = gap_idx

		gap_idx = get_timestamp_gaps_dates(idxss, DIR_out + "/" + t.replace(".pickle", "_gaps_after.html"))

	    savemat(DIR_out+"/summaries_matlab.mat", matlab_output)
	    savemat(DIR_out+"/summaries_matlab_carry_fwd.mat", matlab_output_fwd)

	t1 = time.time()
	print("Time to completion: "+str(datetime.timedelta(seconds=t1-t0)))
