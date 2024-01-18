Pipeline for LAT analysis

Data location

Data can be downloaded at (https://doi.org/10.5281/zenodo.10529801)

Behavioral data and model fits located in Data/Behavior_and_models
Data for testing the fitting procedure located in Data/fitting_pipeline_validation
Neural data located in Data/FR_around_event

Correspondence behavioral and neural data (ID)

Monkey B (Beaker)

181019 = 3
181025 = 8
181029 = 12
181030 = 14
181106 = 18
181109 = 20
181110 = 22
181115 = 26

Monkey S (Scooter)

181102 = 28
181103 = 30
181107 = 32
181108 = 34
181112 = 36
181113 = 38
181116 = 40
181117 = 42
181119 = 44

Reproduce the main figures (all the data necessary is in the corresponding folder)

Figure_1/plot_figure_1/Plot_figure_1_github.m
Figure_2/plot_figure_2/Plot_figure_2_github.m
Figures_3_and_4/plot_figure_3/Plot_figure_3_github_revised.m
Figures_3_and_4/plot_figure_4/Plot_figure_4_github_revised.m
Figure_5/plot_figure_5/Plot_figure_5_github.m
Figure_6/plot_figure_6/Plot_figure_6_github.m

Reproduce the analyses (and the supplementary figures)

You need to set the paths yourself at the start of each script.

Behavioral data (Fig1 - S1-2)

You need to download the VBMC toolbox (https://github.com/lacerbi/vbmc)

In Figure_1/code

doExecuteMasterScriptResetRWVBMC(fsroot,monkey,N_channels,analysis)
save_Model_VBMC
plot_figure1 (behavioral results per monkey. Fig 1 and S1)
plot_models_VBMC (model comparison, Fig S1F)
compare_control_models (model comparison, Fig S1G-H)
plot_models_VBMC_2L (2 learning rate model comparison)

Plot_choice_model_behavior (Fig S2)

To test the fitting pipeline:

Test_parameter_recovery (generate data for each model)
fit_all_models_cluster_proba (for each model and each generated data set)
Test_parameter_recovery_proba (if LOAD=0, will use ‘Fitting_pipeline_validation_results’ in Data/fitting_pipeline_validation to plot Fig S1I-K)

Neural data 

Figure 2 ans S3

In Figure_2/code

Uses:
window [-600 300] around target, size = 900ms
sliding windows around target, -600 to 350ms, size = 300ms
sliding windows around reward end, 0 to 300ms, size = 300ms 

models_FR_figure2_value_function (for window [-600 300])
plot_FigS2_prop_sig_models

prepare_pseudopop_peak_belief_classifier_prog_restricted_FEF
prepare_bootstrap_peak_belief_classifier_across_time_prog
=> same with 4bins and rotated correspond to Fig S3H-I

Run_belief_classifier_combined_time_prog_restricted_FEF (for window [-600 300])
Run_belief_classifier_combined_time_prog_NN_restricted_FEF (for window [-600 300])
=> same with 4bins and rotated correspond to Fig S3H-I, noOne removes neurons with only one block for a template bin

get_projection_classifier_combined_time_prog_restricted_FEF (for window [-600 300])

plot_Fig2_ET_pseudo_pop
=> same with 4bins and rotated correspond to Fig S3H-I, noOne removes neurons with only one block for a template bin

Run_belief_classifier_combined_time_prog_restricted_FEF (for sliding windows)
Run_belief_classifier_combined_time_prog_NN_restricted_FEF (for sliding windows)
get_projection_classifier_combined_time_prog_restricted_FEF (for sliding windows)

plot_Fig2_cross_temporal_decoding

plot_Fig2_PCA

models_FR_figure2_true_template (for window [-600 300])
plot_FigS2_true_template_prop_sig

models_FR_figure2_value_function (for sliding windows)
plot_FigS2_prop_sig_models_across_time

Figures 3, 4 and S4

In Figures_3_and_4/code

Uses:
window [-600 300] around target, size = 900ms

belief_peak_decoder_single_session_update_NN_restricted_FEF.m
Plot_fig3_and_4.m
MDS_template.m

ET_around_reset (Fig S4G)

Figure 5 and S5

In Figure_5/code

Uses:
sliding window around target, -500 to 800ms, size = 200ms
sliding window around response, -500 to 600ms, size = 200ms

prepare_pseudopop_choice_prog_classifier_across_time
prepare_bootstrap_peak_belief_choice_across_time
Run_classifier_peak_belief_choice_restricted_FEF
Accuracy_classifier_peak_belief_choice

prepare_pseudopop_peak_belief_value_classifier_across_time
prepare_bootstrap_peak_belief_value_across_time
Run_classifier_peak_belief_value_restricted_FEF
Accuracy_classifier_peak_belief_value

prepare_pseudopop_peak_belief_cc_classifier_across_time
prepare_bootstrap_peak_belief_cc_across_time
Accuracy_classifier_peak_belief_chosen_color_2bins

prepare_pseudopop_stim_color_across_time
prepare_bootstrap_stim_color_across_time
Run_stim_color_classifier_across_time_with_NN_restricted_FEF
Accuracy_classifier_stim_color

Figure 6 and S6

In Figure_6/code

Uses:
sliding window around target, -400 to 800ms, size = 200ms

get_explained_variance_split_restricted_FEF
get_boot_glm_all_values_split_correct
get_mean_EV_across_split_and_locs
get_preferred_loc_restricted_FEF
get_boot_glm_all_values_split_pref

plot_correlation_chosen_unchosen
plot_correlation_local_global
plot_correlation_across_locations
plot_EV_across_loc
plot_correlation_chosen_unchosen_pref







