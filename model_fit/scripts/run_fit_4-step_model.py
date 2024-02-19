# script to run the fit for the 4-step model including nuclear decay

# --- PACKAGES ----------------------
from pathlib import Path
import numpy as np
import pandas as pd
from lmfit import Parameters, Model, minimize
rng = np.random.default_rng()
# import model functions from helper file
from general_solutions_4_step_model import *
# import the fit parameters
from fit_parameters import *


# --- select gene to fit from snakemake input / wildcard ----------
gene_id = snakemake.wildcards.gene



##############################################################################
# --- DATA IMPORTS --------------------------------------------------
##############################################################################

# --- steady state data ---------
ss_gene = (pd.read_table(snakemake.input.steadystate)
 .filter(['cyto_over_nuc', 'se_cyto_over_nuc', 
          'mem_over_cyto', 'se_mem_over_cyto'])
 .iloc[0]
 .to_numpy()
)

# --- T2C data ------------------
# import T2C data; exons with low T counts are filtered out;
# compartments are ordered as: nuc_exon, cyto, mem, nuc_intron
# data is ordered by exons (gene-level intron data comes last)
T2C_data = pd.read_table(snakemake.input.T2C)

# --- exon distance to end and gene length ------------
# import start, stop and distance to 3' end of each exon
start_stop = pd.read_table(snakemake.input.start_stop)

# --- gene length in kb -----------
gene_length = abs(start_stop.stop.max() - start_stop.start.min()) / 1000
# --> as time values are in min, this makes the elongatoin rate kb/min

# --- dataframe with T2C and distance to end data --------------------
# merge T2C and start_stop data
all_data = T2C_data.merge(start_stop.filter(['exon_id', 'distance_to_end_in_kb']), 
                      how='left')
# add gene length in distance to end column
all_data.loc[:,'distance_or_length'] = np.where(
    all_data['exon_id'] == 'gene',
    gene_length,
    all_data['distance_to_end_in_kb']
)



##############################################################################
### MODELING
##############################################################################

# --- FUNCTIONS TO EVALUATE MODEL -------------------------
# we have the functions/solutions of compartments imported
# now make functions that calculate the fit dataset for each compartment
def nuc_pre_dataset(params, t, delay, gene_length):
    k1 = np.exp(params['k1'])
    elongation_rate = np.exp(params['elongation_rate'])
    t_mod = modified_time(t, elongation_rate, gene_length/10, delay)
    return step1(t_mod, k1)

def nuc_mat_dataset(params, t, delay, distance_to_end):
    k1 = np.exp(params['k1'])
    k2 = np.exp(params['k2'])
    g2 = np.exp(params['g2'])
    elongation_rate = np.exp(params['elongation_rate'])
    t_mod = modified_time(t, elongation_rate, distance_to_end, delay)
    return step2(t_mod, k1, k2+g2)

def cyto_dataset(params, t, delay, distance_to_end):
    k1 = np.exp(params['k1'])
    k2 = np.exp(params['k2'])
    g2 = np.exp(params['g2'])
    k3 = np.exp(params['k3'])
    # g3 = np.exp(params['g3'])
    elongation_rate = np.exp(params['elongation_rate'])
    t_mod = modified_time(t, elongation_rate, distance_to_end, delay)
    return step3(t_mod, k1, k2+g2, k3)

def mem_dataset(params, t, delay, distance_to_end):
    k1 = np.exp(params['k1'])
    k2 = np.exp(params['k2'])
    g2 = np.exp(params['g2'])
    k3 = np.exp(params['k3'])
    # g3 = np.exp(params['g3'])
    g4 = np.exp(params['g4'])
    elongation_rate = np.exp(params['elongation_rate'])
    t_mod = modified_time(t, elongation_rate, distance_to_end, delay)
    return step4(t_mod, k1, k2+g2, k3, g4)

#### define function that calculates residual from dataframe
def calc_residuals_from_df(params, df, delay):
    t_array = df['time'].to_numpy()
    dist_array = df['distance_or_length'].to_numpy()
    # evaluate fcuntions with parameters
    df.loc[:,'function_eval'] = pd.to_numeric(np.select(
        [
            df['fraction'] == 'nuc_intron',
            df['fraction'] == 'nuc_exon',
            df['fraction'] == 'cyto',
            df['fraction'] == 'mem',
        ],
        [
            nuc_pre_dataset(params, t_array, delay, dist_array),
            nuc_mat_dataset(params, t_array, delay, dist_array),
            cyto_dataset(params, t_array, delay, dist_array),
            mem_dataset(params, t_array, delay, dist_array),
        ],
        np.nan
    ))
    # calcualte residuals with weights
    df = (df
     .assign(residuals_t2c = lambda x: (x['T2C_over_T'] - x['function_eval']) / x['weighted_error'])
    )
    res = df['residuals_t2c'].to_numpy()
    return res

# define function that calculates total flattened residual
def objective(params, df, delay):
    resid = np.zeros(len(df) + 2)
    # --- make residuals per compartment -------------
    resid[:len(df)] = calc_residuals_from_df(params, df, delay)

    # add residuals for steady state
    resid[-2] = (ss_gene[0] - params['cyto_over_nuc']) / ss_gene[1]
    resid[-1] = (ss_gene[2] - params['mem_over_cyto']) / ss_gene[3]
    return resid



##############################################################################
### FIT LOOP
##############################################################################

# create empty lists to append the result of each loop to
list_values = []
list_chis = []

# loop to run multiple times with random initial values
for i in range(samples):
    # ---- CREATE FIT PARAMETER OBJECT ----------------------------
    # the nice way (use errors of steady state)
    # create parameter object
    pars = Parameters()
    # add normal parameters
    pars.add('k1', value=random_parameter(lower_start, upper_start),
             min=lower_bound, 
             max=upper_bound)
    pars.add('k2', value=random_parameter(lower_start, upper_start),
             min=lower_bound, max=upper_bound)
    pars.add('g2', value=random_parameter(lower_start_decay, upper_start),
              min=lower_bound_decay, 
              max=upper_bound)
    pars.add('k3', value=random_parameter(lower_start, upper_start),
             min=lower_bound, 
             max=upper_bound)
    # pars.add('g3', value=random_parameter(lower_start, upper_start),
    #          min=lower_bound_decay, 
    #           max=upper_start)
    pars.add('g4', value=random_parameter(lower_start, upper_start),
             min=lower_bound, 
             max=upper_bound)
    # constrain them using expr and min and max
    pars.add('cyto_over_nuc', expr='exp(k2) / exp(k3)',
             min=ss_gene[0]-5*ss_gene[1], max=ss_gene[0]+5*ss_gene[1])
    pars.add('mem_over_cyto', expr='exp(k3) / exp(g4)',
             min=ss_gene[2]-5*ss_gene[3], max=ss_gene[2]+5*ss_gene[3])
    # add elongation rate
    pars.add('elongation_rate', value=random_parameter(ls_erate, us_erate),
             min=lb_erate, max=ub_erate)

    # --- RUN THE FIT --------------------
    # actual fitting
    out = minimize(objective, pars, args=(all_data, overall_delay), 
                   method=method, scale_covar=False, calc_covar=False, 
                   nan_policy='omit')

    # get fitted parameters
    k1 = np.exp(out.params['k1'].value)
    k2 = np.exp(out.params['k2'].value)
    g2 = np.exp(out.params['g2'].value)
    k3 = np.exp(out.params['k3'].value)
    g3 = np.nan
    g4 = np.exp(out.params['g4'].value)
    v = np.exp(out.params['elongation_rate'].value)

    # --- make values that best determine fit quality (besides chi-squared) ------
    # calculate how close fit values are to upper boundaries ("boundary cost")
    # if they are too close, they are unbiological --> filter them out later
    # only export or decay rates
    b_cost_export = (
                    1/(np.exp(upper_bound) - k1 + offset)**2
                    + 1/(np.exp(upper_bound) - k2 + offset)**2
                    + 1/(np.exp(upper_bound) - g2 + offset)**2
                    + 1/(np.exp(upper_bound) - k3 + offset)**2
                    # + 1/(np.exp(upper_bound) - g3 + offset)**2
                    + 1/(np.exp(upper_bound) - g4 + offset)**2) / 16
    # elongation rates with difference to upper and lower boundary
    b_cost_elong = (1/(np.exp(ub_erate) - v + offset)**2 / 10
                   + 1/(v - np.exp(lb_erate) + offset)**2 / 16)

    # create list of result from fit
    values = [k1, k2, g2, k3, g3, g4, v]
    chis = [out.redchi, b_cost_export, b_cost_elong]

    # --- append the results of each fit ---------------
    list_values.append(values)
    list_chis.append(chis)



##############################################################################
### GET BEST FIT OUTPUT VALUES AND SAVE ALL RESULTS
##############################################################################

# ---- create dataframe with all values and chis ------
df_values = pd.DataFrame(np.array(list_values))
df_values.rename(columns={0: 'nuc_intron_export', 1: 'nuc_exon_export',
                          2: 'nuc_exon_decay', 
                          3: 'cyto_exon_export', 4: 'cyto_exon_decay',
                          5: 'mem_exon_decay', 6: 'elongation_rate', },
                          inplace=True)

df_chis = pd.DataFrame(np.array(list_chis))
df_chis.rename(columns={0: 'red_chi',
                        1: 'b_cost_export', 2: 'b_cost_elong'}, inplace=True)

all_results = df_values.join(df_chis)
all_results = all_results.sort_values('red_chi')
# write out results from all samples to gene's subdirectory
all_results_path = Path(str(snakemake.output))
all_results.to_csv(all_results_path, sep='\t', index=False)


# --- get the best fit result -----
# make df with top ten fit results
top_ten_results = all_results.iloc[:10,:]
# get smallest red chi-squared value and ad small tolerance
cutoff_red_chi = top_ten_results['red_chi'].min() + 0.5
# get index of best result based on small red chi with smallest b_cost
best_fit_index = (top_ten_results
 .query('red_chi < @cutoff_red_chi')['b_cost_export']
 .idxmin()
)

# get values
best_fit = all_results[all_results.index == best_fit_index]    # returns dataframe
best_fit.insert(0, 'gene_id', gene_id)

# append best fit result to final table
output_path = Path('fit_table_until_membrane_no_cyto_decay.tsv')
best_fit.to_csv(output_path, mode='a', sep='\t',
                header=not output_path.exists(), index=False)


# --- errors from chi sqr statisitcs ---------
# take top ten fit results -> calculate standard deviation of good fit results
top_ten_results.insert(
    3, 'nuc_exon_removal', 
    top_ten_results['nuc_exon_export'] + top_ten_results['nuc_exon_decay']
)
errors = top_ten_results.std().to_frame().T
errors = errors.loc[:, 'nuc_intron_export':'elongation_rate']
errors.columns = ['std_' + x for x in errors.columns]
errors.insert(0, 'gene_id', gene_id)
# append deviation / error result to final error table
error_path = Path('standard_deviation_until_membrane_no_cyto_decay.tsv')
errors.to_csv(error_path, mode='a', sep='\t',
              header=not error_path.exists(), index=False)
