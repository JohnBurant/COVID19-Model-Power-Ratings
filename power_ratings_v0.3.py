# power_ratings_v0.2.py
# Last updated: 22 August 2020
# https://github.com/JohnBurant/COVID19-Model-Power-Ratings

# Timeframe choices
first_week = 17
last_week = 33
nweeks = 3 # Input as 0-indexed, i.e., how many we actually want - 1
Nroll = 4 # Number of weeks to use for rolling average overall power rating (N->inf = lifetime, N->1 identical to weekly rating)

# Policy choices
# See discussion on github site
d_policy = {'include_baselines': True,
            'include_ensemble' : True,
            'include_models_without_full_measure_set': False,
            'include_models_without_all_weeks': False,
            'generate_partial_results_for_recent_weeks': False
           }

working_dir = 'EDIT_ME_PLEASE' # Local path of this repo
yyg_repo_dir = 'EDIT_ME_PLEASE' # path to clone of https://github.com/youyanggu/covid19-forecast-hub-evaluation
eval_dir = yyg_repo_dir + 'evaluations/'
output_dir = working_dir + 'results/'

import pandas as pd
import numpy as np
import epiweeks
import datetime
import os
import altair as alt
idx = pd.IndexSlice

def power_score(s_in):
    '''
    Given a pd.Series of error measures, return a (sorted) pd.Series w/ power ratings.
    Power rating is as follows:
      The minimum error is assigned a power rating of 100;
      The median error is assigned a power rating of 50;
      Remaining errors are distributed proportionately according to how they fall on the scale defined by min->100, median->50;
      Any errors that would result in a power rating < 0 are set = 0.
    '''
    s = s_in.sort_values()
    err_range = 2*(s.median()-s.min())
    s_scaled = (err_range - s + s.min()).clip(0)
    s_scaled = 100*s_scaled/s_scaled[0]
    return s_scaled


def say(st,l_outs,save=True,term='\n'):
    '''
    Given a string st, a list l_outs, a boolean save and an options string term:
      print st, using term asthe argument end of print function
      if save=True append st to l_outs
    '''
    print(st,end=term)
    if save: l_outs.append(st)

l_outs = []

say('Beginning generation of power ratings from raw evaluation files.',l_outs)
datestamp = datetime.datetime.today().strftime('%Y-%m-%d')
outstring = 'Datestamp: '+datestamp
say(outstring,l_outs)
say('Following these policies:',l_outs)
for k, v in d_policy.items():
    outstring = '  '+str(k)+': '+str(v)
    say(outstring,l_outs)


# Create some handy dicts to move between YYYY-MM-DD and epiweek
d_epiweeks = {i: epiweeks.Week(2020,i) for i in range(first_week,last_week+1)}
d_wk_to_enddate = {k: w.enddate().strftime('%Y-%m-%d') for k, w in d_epiweeks.items()}
d_enddate_to_wk = {v: k for k, v in d_wk_to_enddate.items()}
d_wk_to_startdate = {k: (w.startdate()+datetime.timedelta(days=1)).strftime('%Y-%m-%d') for k, w in d_epiweeks.items()}
d_startdate_to_wk = {v: k for k, v in d_wk_to_startdate.items()}

# For each starting week, create list of week pairs that are (starting_week,ending_week). Store as dict.
d_week_pairs = {i: [(j,i) for j in range(first_week,last_week+1) if i >= j] for i in range(first_week,last_week+1)}

# Also do this with the full string YYYY-MM-DD version
d_week_pairs_str = {}
for k, v in d_week_pairs.items():
    d_week_pairs_str[d_wk_to_enddate[k]] = [d_wk_to_startdate[tup[0]]+'_'+d_wk_to_enddate[tup[1]] for tup in v]

# Read in all the files in the evaluations directory
# We could just glob but curating it slightly might make the manipulations easier
# Note we make the switch from YYYY-MM-DD to integer weeks right here.

# Format {filename_base: {'cols': list of columns to keep, 'short_name': what we call the measure}
d_files_measures_cols = {'states_abs_errs': {'cols' : ['mean'], 'short_name': 'state_abs'},
                         'states_sq_errs': {'cols' : ['mean'], 'short_name': 'state_sq'},
                         'states_mean_ranks': {'cols' : ['mean_rank'], 'short_name': 'state_rank'},
                         'us_errs': {'cols': ['perc_error'], 'short_name': 'us_tot'}
                        }

measures_used = d_files_measures_cols.keys()
measures_used_str = 'Measures used: '+', '.join(list(measures_used))
say(measures_used_str,l_outs)
say('Loading evaluation files...',l_outs,term='')
l_dfs = []
l_files_not_found = [] # Can inspect this later if we're curious
ict = 0
for meas, vd in d_files_measures_cols.items():
    cols_to_keep = vd['cols']
    short_name = vd['short_name']
    for k, v in d_week_pairs_str.items():
        for pr in v:
            ict += 1
            fstring = k+'/'+pr+'_'+meas+'.csv'
            full_path_to_file = eval_dir+fstring
            # Before we try to read the file check if it exists?
            if os.path.isfile(full_path_to_file):
                df = pd.read_csv(full_path_to_file)
                df = df[['Unnamed: 0', *cols_to_keep]]
                df = df.set_index('Unnamed: 0')
                date_pr = pr.split('_')
                int_weeks = [d_startdate_to_wk[date_pr[0]], d_enddate_to_wk[date_pr[1]]]
                cols_mi = pd.MultiIndex.from_tuples([(short_name+'_'+c, int_weeks[0], int_weeks[1]) for c in cols_to_keep],
                                               names=['measure','start_wk','end_wk'])
                df.columns = cols_mi
                l_dfs.append(df)
            else:
                l_files_not_found.append(full_path_to_file)

say('Done',l_outs,save=False)
outstring = 'Number of possible week-pair files: '+str(ict)
say(outstring,l_outs)
outstring = 'Number of files processed: '+str(len(l_dfs))
say(outstring,l_outs)

df_raw = pd.concat(l_dfs,axis=1)
df_raw = df_raw.dropna(axis=0,how='all') # Some models show up but not for the metric we use (notably: UChicago, which did only Illinois)
df_raw.index.name = 'Model'
outstring = 'Shape of df_raw: '+str(df_raw.shape)
say(outstring,l_outs)

# Do some further processing into a new df
#   df will be the main dataframe we'll work with that will contain start_wk/num_wk/model/measure (4-tuple) error measures
# Switch from first_wk/end_wk to first_wk/num_wks where num_wks counts how many weeks out from first_wk the predictions are, starting w/ 0.
df = df_raw.copy()
df = df.T.reset_index()
df['num_wks'] = df['end_wk'] - df['start_wk']
df = df.drop(columns=['end_wk'])
df = df.set_index(['start_wk','num_wks','measure'])
df = df.sort_index()

# Also, slight processing:
# (1) convert mean-square-errors into root-mean-square-errors (and rename the index)
# (2) take absolute value of percentage error on whole-US projection

# (1a) apply the sqrt
screen_sq_errs = df.index.get_level_values('measure').str.startswith('state_sq')
df.loc[screen_sq_errs,:] = df.loc[screen_sq_errs,:].applymap(np.sqrt)

# (1b) rename the measure; easier to do this way than dealing with the index object itself!
df = df.reset_index('measure')
df.loc[screen_sq_errs,'measure'] = 'states_rmse'
df = df.set_index('measure',append=True)

# (2)
screen_perc_errs = df.index.get_level_values('measure') == 'us_tot_perc_error'
df.loc[screen_perc_errs,:] =  df.loc[screen_perc_errs,:].applymap(lambda x: x if pd.isna(x) else abs(float(x.strip('%'))))

# Also, a slightly nicer name for the state rank
df = df.reset_index('measure')
screen_state_rank = df['measure'] == 'state_rank_mean_rank'
df.loc[screen_state_rank,'measure'] = 'state_mean_rank'
df = df.set_index('measure',append=True)

# Now that we've fixed percentages, make everything a float
df = df.astype(float)

# Trim the weeks further into the future from the start_wk than we want to use, and sort the index for efficiency later
screen = df.index.get_level_values('num_wks') <= nweeks
df = df.loc[screen,:].sort_index()

# Optionally, take out baseline(s) COVIDhub ensemble from the power rating process depending on policy choices
l_excluded_models = ['Baseline'] # Always exclude YYG's 'Baseline' b/c we already have another baseline (COVIDhub-baseline) in here
if not d_policy['include_baselines']: l_excluded_models += ['COVIDhub-baseline']
if not d_policy['include_ensemble']: l_excluded_models += ['COVIDhub-ensemble']
cols = [c for c in df.columns if c not in l_excluded_models]
df = df[cols]

df.columns.name = 'Model' # This isn't carried through manipulations above?

outstring = 'Generating working df of shape: '+str(df.shape)
say(outstring,l_outs)

# Calculate the power ratings for each start_wk/num_wks/measure/Model 4-tuple, where the power rating is done over all models present for the start_wk/num_wks/measure 3-tuple
# df_power holds the power ratings
df_power = df.apply(power_score,axis=1)
df_power.columns.name = 'Model' # power_score doesn't know about name of series it returns

# Make some aggregations depending on policy choices
# dft is modified and used throughout; we make some worthwhile intermediate aggregations and store as we go.

# First, create a single power rating for each Model/start_wk/num_wks tuple (averaging over all measures)
dft = df_power.stack().unstack('measure')
if d_policy['include_models_without_full_measure_set']:
    dft = dft.fillna(0,axis=1)
dft = dft.mean(skipna=False,axis=1) # It's actually a series at this point
df_power_model_wk_pair = dft.unstack('Model')

# Now, for each Model/start_wk 2-tuple, aggregate over the weeks (up to nweeks)

# First, count the number of weeks including/after start_wk for which there is any data for each model
s_model_wks = dft.unstack(level=['num_wks']).apply(lambda x: x.notna().sum(),axis=1)

# Now for each start week, how many weeks data are expected?
i = s_model_wks.index.get_level_values('start_wk').unique()
max_weeks_expected = nweeks + 1 - np.maximum(i-last_week+nweeks,0)
s_start_wk_max_wks = pd.Series(index=i,data=max_weeks_expected)

# Now check if the weeks we counted as having data is sufficient
screen_sufficient_weeks = (s_model_wks.unstack('Model').fillna(0).astype(int)
                           .apply(lambda x: x.values == s_start_wk_max_wks.values).stack('Model'))
i_model_wks = s_model_wks[screen_sufficient_weeks].index

# Finally, apply the policies and aggregation
dft = dft.unstack('num_wks')
if not d_policy['include_models_without_all_weeks']:
    dft = dft.loc[i_model_wks,:]
dft = dft.mean(skipna=d_policy['generate_partial_results_for_recent_weeks'],axis=1).dropna(axis=0,how='all')

df_power_model_start_wk = dft.unstack('Model')
df_power_model_start_wk.index.name = 'start_wk'
df_power_model_start_wk.columns.name = 'Model'

# Create the rolling mean, and change index to be a little more intuitive
df_power_model_rolling_N_wks = df_power_model_start_wk.rolling(Nroll).mean().dropna(axis=0,how='all').T
df_power_model_rolling_N_wks.columns = [str(c-4+1)+'-'+str(c) for c in df_power_model_rolling_N_wks.columns]
df_power_model_rolling_N_wks.columns.name = 'start_wk_range'

outstring = 'Using '+str(Nroll)+' weeks window for mean rolling power rating'
say(outstring,l_outs)

# Generate a Model's lifetime average power ranking, over a specified set of weeks (can be all weeks of course)
week_alpha = first_week
week_omega = max(df_power_model_start_wk.index)+1
weeks_included = list(range(week_alpha,week_omega))
outstring = 'Generating an average of power ratings for models over a span of '+str(len(weeks_included))+' weeks.'
say(outstring,l_outs)
outstring = 'First week: '+str(first_week)+' (i.e., week starting '+str(d_wk_to_startdate[first_week])+')'
say(outstring,l_outs)
outstring = '  Covering the period of first week plus an additional '+str(nweeks)+' weeks (i.e., a total of '+str(nweeks+1)+' weeks)'
say(outstring,l_outs)
df_power_model = df_power_model_start_wk.loc[weeks_included,:].mean(axis=0).sort_values(ascending=False)
df_power_model.name = 'Power_Rating'

# Also look at national-only models
has_national = df_power.loc[idx[:,:,'us_tot_perc_error'],:].notna().droplevel('measure')
has_all_state = df_power.loc[idx[:,:,['state_abs_mean','state_mean_rank','states_rmse']],:].notna().groupby(['start_wk','num_wks']).sum().astype(int)==3
is_national_only = has_national & ~ has_all_state
national_only_week_has_4plus_weeks = (is_national_only.groupby('start_wk').sum().astype(int)>0).T.sum(axis=1).astype(int)>4
l_national_only_models = national_only_week_has_4plus_weeks[national_only_week_has_4plus_weeks].index.to_list()

# Build the power measure df for just the national measures; filter according to policy and add indication about whether a model is national-only
screen_us_meas = df_power.index.get_level_values('measure')=='us_tot_perc_error'
dfx = df_power[screen_us_meas].droplevel('measure') # .groupby('start_wk') .mean().T

# Note: This code is duplicated above when used for the same purpose to build df_power_model_start_wk
s_model_wks = dfx.stack().unstack('num_wks').apply(lambda x: x.notna().sum(),axis=1)

# Now for each start week, how many weeks data are expected?
i = s_model_wks.index.get_level_values('start_wk').unique()
max_weeks_expected = nweeks + 1 - np.maximum(i-last_week+nweeks,0)
s_start_wk_max_wks = pd.Series(index=i,data=max_weeks_expected)

# Now check if the weeks we counted as having data is sufficient
screen_sufficient_weeks = (s_model_wks.unstack('Model').fillna(0).astype(int)
                           .apply(lambda x: x.values == s_start_wk_max_wks.values).stack('Model'))
i_model_wks = s_model_wks[screen_sufficient_weeks].index

dfx = dfx.stack().unstack('num_wks')
if not d_policy['include_models_without_all_weeks']:
    dfx = dfx.loc[i_model_wks,:]
dfx = dfx.mean(skipna=d_policy['generate_partial_results_for_recent_weeks'],axis=1).dropna(axis=0,how='all')

df_power_national_measure = dfx.unstack('Model')
df_power_national_measure.index.name = 'start_wk'
df_power_national_measure.columns.name = 'Model'
df_power_national_measure = df_power_national_measure.T

df_power_national_measure['national_only_model'] = df_power_national_measure.index.isin(l_national_only_models)
df_power_national_measure = df_power_national_measure.set_index('national_only_model',append=True)

# We can also look at models over their whole lifetime (or, of course, a segment of it) and calculate each model's power rating for each specific measure.
df_power_model_measure = df_power.stack().groupby(['Model','measure']).mean().unstack('measure')

# We can also look at the range (e.g., max-min) of the power ratings of a model for the various measures
df_power_model_measure.dropna(axis=0,how='any').apply(lambda x: x.max() - x.min(),axis=1).sort_values()

# We could also look at power rating by number of weeks out from initial forecast (i.e., agg by num_wks)
df_power_model_num_wks = df_power.stack().groupby(['Model','num_wks']).mean().unstack()

# Save the work.
output_filename_base = "covid_model_power_ratings"
output_file_week_str = 'start_wk_'+str(first_week)+'-num_wks_'+str(nweeks+1)
output_filename = output_filename_base+'-'+output_file_week_str+'-'+datestamp+'.xlsx'
full_path_to_output_file = output_dir+output_filename

# (Possibly more) complicated (than necessary) way to set up the Info worksheet
d_output_sheets = {1: {'df': df_power_model,
     'sheet_name': 'Model_lifetime',
     'desc': 'Power rating for each model over its full lifetime (single number for each model)'},
 2: {'df': df_power_model_start_wk.T, # the transpose makes for better viewing
     'sheet_name': 'Model_weekly',
     'desc': 'Power ratings for each model for each week projections were submitted (single number for each model for each week)'},
 3: {'df': df_power_model_rolling_N_wks,
     'sheet_name': 'Model_rolling_N_wks',
     'desc': 'Power ratings for each model averaged over a rolling N week window, N given in Info tab'},
 4: {'df': df_power_model_measure,
     'sheet_name': 'Model_measures',
     'desc': 'Power ratings for each separate measure for each model over its full lifetime (n_measures number for each model)'},
 5: {'df': df_power_model_num_wks,
     'sheet_name': 'Model_num_wks',
     'desc': 'Power ratings for each model for increasing 0, 1, 2, etc. (0-indexed) week periods since a projection was submitted, averaged over that projection\'s full lifetime (n_weeks numbers per model)'},
 6: {'df': df_power_model_wk_pair,
     'sheet_name': 'Model_week_pairs',
     'desc': 'Power ratings for each model for increasing 0, 1, 2, etc. (0-indexed) week periods since a projection was submitted, for each week projections were submitted'},
 7: {'df': df_power_national_measure,
     'sheet_name': 'Model_natl_err_only_wks',
     'desc': 'Power rating for the national death total error measure only; indicates if a model is a multi-state or national-only model'},
 8: {'df': df_power,
     'sheet_name': 'Raw_power_ratings',
     'desc': 'Power ratings for every model, for each model/start_wk/num_wk/measure 4-tuple'}
}
l_sheet_names = [vd['sheet_name'] for vd in d_output_sheets.values()]
l_descs = [vd['desc'] for vd in d_output_sheets.values()]

df_index_sheet = pd.DataFrame(data=[l_sheet_names,l_descs],index=['Sheet name','Description']).T

outstring = 'Saving power ratings in '+output_filename+'...'
say(outstring,l_outs,term='')

df_info_sheet = pd.DataFrame(l_outs)

with pd.ExcelWriter(full_path_to_output_file) as writer:
    df_info_sheet.to_excel(writer,sheet_name='Info',header=False,index=False)
    df_index_sheet.to_excel(writer,sheet_name='Index',header=True,index=False)
    #df_info.to_excel(writer, sheet_name='Info',header=False,index=False)
    for k, v in d_output_sheets.items():
        sheet_name = v['sheet_name']
        df = v['df']
        df.to_excel(writer, sheet_name=sheet_name,float_format='%.2f')
say('Done',l_outs,save=False)
