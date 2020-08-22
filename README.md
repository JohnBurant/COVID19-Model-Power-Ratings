# COVID19 Model Power Ratings

COVID19 Model Power Ratings provide a framework for comparing the accuracy of the models submitted to @reichlab's COVID19 Forecast Hub (https://covid19forecasthub.org/) as a ***single number***.

At present, the power ratings are built from the error measures calculated by @youyanggu's Evaluation of COVID-19 Models (https://github.com/youyanggu/covid19-forecast-hub-evaluation/) package. That package provides a variety of error measures, but does not consolidate them into a single number. The measures in @youyanggu's package measure only the accuracy of forecasted deaths and therefore the power ratings here are also limited to forecasted deaths. A power ratings calculations does not, by itself, require the calculation of any new error measures. In fact the code presented here simply reads the error measures contained in @youyanggu's covid19-forecast-hub-evaluation output files and uses them as input to further calculations.

The power rating framework provides the flexibility to compare models over any desired time range and time horizon, to include or exclude any particular models, and to use  any desired set of error measures. This flexibility enables both high-level summary ratings and more detailed analyses of the accuracy of the models.

### Desired Attributes of a Comparative Measure ###

A method that compares the accuracy of these models over time should:
* Provide a consistent numerical scale from week to week.
* Provide an output numerical measure that describes only the models' performance relative to each other, not to absolute performance. (Although models are assessed on the basis of their absolute performance as input the calculation of comparative measures, the final numbers returned as comparative measures should depend only *indirectly* on the absolute performance. Stated another way: If all models perform equally relatively well or relatively poorly on an absolute measure, the measure of their relative performance should not change.)
* Consider performance over the full time range of a forecast's time horizon, e.g., if the time horizon of a forecast is 4 weeks, the evaluation should consider not just the 4-week period, but all intermediate periods (1-week, 2-week, 3-week and 4-week) periods over that time horizon.
* Consider multiple geographical levels of aggregation, e.g., both national & state level. It should take into account that average error comparisons across multiple geographies are by definition skewed towards performance in the most populous geographies (technically those geographies with the highest number of deaths).
* Distinguish between relatively close performance and relatively disparate performance (i.e., it should not be a simple ranking).

### Power Rating ###

The power rating framework has these attributes.

For any particular *individual* error measure applied to a set of models, the power rating scoring is as follows:

* The model with the minimum error is assigned a power rating of 100.
* The model with the median error is assigned a power rating of 50.
* Other models are assigned a power rating on a proportional basis, e.g.
  * A model with an error midway between the minimum error and the median error is assigned a power rating of 75.
  * A model with an error 1.5x the difference between the minimum error and the median error is assigned a power rating of 25.
  * A model with an error 2x the difference between the minimum error and the median error is assigned a power rating of 0.
* A power rating score cannot be less than 0. Therefore any model with an error >2x the difference between the minimum error and the median error is also assigned a power rating of 0.

Truncating the low-end ratings to be 0 at 2x the difference between the minimum error and the median error ensures that the power ratings are not skewed by the inclusion of highly inaccurate (outlier) models in the set.

### Choices ###

The power rating as described above is a generic scheme that could be used in many settings. Applying it to the models on the COVID19 Forecast Hub site, however, requires making choices of certain parameters and inclusion/exclusion policies.

*Parameters*

Time frame: Model forecasts included in the COVID19 Forecast Hub are submitted and evaluated on a weekly cycle. The earliest forecasts date to week 17 of 2020 (the week ending 2020-04-25). By default, the code here includes all forecasts from this week forward.

Forecast time horizon: Forecasts submitted to the COVID19 Forecast Hub typically have a time horizon of at least four weeks from the date of submission (including forecasts for that same week plus three additional weeks). Some models includes forecasts for a longer time horizon. The default time horizon used by the code here is 4 weeks.

*Measures Included*

The default set of measures taken from @youyanggu's Evaluation of COVID-19 Models package are:

* Mean absolute error (state-level)
* Root-mean-square error (state-level) (note: @youyanggu provides square errors, the code here takes the square root of these)
* Mean ranking (state-level)
* Absolute percent error (national)

By virtue of including three state-level measures, the default power rating emphasizes state-level accuracy over national-level accuracy. Because of the disparate levels of impact of COVID19 across the country and the value of having localized forecasts, this is a reasonable choice. The highly uneven population distribution among the states, however, means that any average of state-level measures will be skewed towards performance in the more populous states.
The inclusion of the mean ranking (i.e., the averaged ranking of the models in each of the states) rewards consistency across states and offsets this issue somewhat.
The inclusion of the absolute error at the national level provides credit for models that may not be highly accurate for all states, but for which cancellation of errors leads to accuracy at the national level.

*Model Inclusion/Exclusion*

The power ratings are inherently dependent on the set of models included in the ratings. Several options are permitted by the current code:

* **UPDATED** Whether or not to include baseline models (linear projections). The default is to include **one baseline model**.
* Whether or not to include the COVID19 Forecast Hub ensemble model, which itself is an aggregation of the other models. The default is to include the ensemble.
* **UPDATED** Whether or not to include models for which the full set of error measures are not available. Most models provide projections at the state level, which are then aggregated to the national level. However, some models only provide national forecasts, and thus are missing 3 of the 4 error measures. Models that provide only national forecasts receive low overall power scores if compared with state-level models (as they are assigned a 0 for missing measures), so it is reasonable to analyze the national-only models separately and not to include them on lists of state-level model power ratings. However, it is still worthwhile to compare the aggregation of state-level models against the full set of state-level and national-only models. **This is now the default setting.**
* Whether or not to include models that do not provide forecasts for all weeks in the desired time horizon. To compare models with each other consistently, it's reasonable to require models to have projections available for the identical number of weeks. Therefore the default is not to include models that do not provide projections for all weeks in the desired time horizon.
* Whether or not to generate results for recent weeks. The power ratings are assigned on the basis of performance over a time horizon of specified number of weeks. This policy choice determines whether (preliminary) power ratings are generated for model forecasts prior to reaching the end of that time horizon. The default is not to generate these preliminary ratings.

### Process ###

The first power ratings calculated are for each start_week/num_week/measure/model 4-tuple (where num_week is the 0, 1, 2, etc. week in the desired time horizon). To calculate these power ratings, all models are compared across every start_week/num_week/measure 3-tuple.

The 4-tuple power ratings are then aggregated by straight averaging in a variety of ways to calculate the aggregated power ratings.

The feature measures of this code are:

**Lifetime power rating**: Provides a single number to compare each model with all other models across its whole lifetime, for the chosen set of error measures and time horizon.

**Weekly power rating**: Provides a weekly measure for each model, for each week in which a model forecast was provided, comparing it against all other models for the chosen set of error measures and time horizon.

**UPDATED** **Rolling average power rating**: Provides a measure of which models have been accurate at specific points in the past and which are the most accurate *in the current moment*. This calculates a rolling average of the weekly power ratings over a window (default: four week window). The four-week averaging reduces the flucuation realtive to the purely weekly power rating, and therefore allows more meaningful assessment over periods of a month, as well as demonstrating more clearly overall trend lines.

### Code ###

The power rating scheme is implemented in Python. The code can be found in **(UPDATED) power_ratings_v0.3.py**, in the root directory of this repo. The code requires a local copy of the /evaluation subdirectory of https://github.com/youyanggu/covid19-forecast-hub-evaluation.

### Results ###

The python code generates a multi-tab Excel file, which can be found in the /results subdirectory of this repo.

*Info:* Summarizes parameters used in the power rating calculation and provides a datestamp.

*Model_lifetime:* The lifetime power ratings as described above.

*Model_weekly:* The weekly power ratings as described above.

**UPDATED** *Model_rolling_N_wks:* Provides the rolling average power rating over N weeks (N listed in info tab) as described above.

*Model_measures:* Provides a comparison of each model on each error measure over its whole lifetime. This shows which models do well, or not well, on each measure. Some are relatively consistent across all measures, others vary widely.

*Model_num_weeks:* Provides a comparison of each model by the number of weeks out from the week of submission. This shows how well models perform over the course of the time horizon. Some models improve over the course of the time horizon, whereas others get worse.

*Model_week_pairs:* Provide a power rating for a model for each week throughout the time horizon, separately for each week in which a forecast was submitted.

**UPDATED** *Model_natl_err_only_wks:* Power ratings of all models on the basis of the error in national death projections only, provided separately for each week a forecast was submitted.

*Raw_power_ratings:* The power ratings for the 4-tuple of model/start_week/end_week/measure described above, which form the basis for all power rating aggregations.
