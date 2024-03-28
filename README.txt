In the following folders, we conduct simulations with sigma = 100 and psi = 0.01:

	Fam-simulation-200/Miss-20: Simulation with I = 200 and censoring rate of 20%
	Fam-simulation-200/Miss-50: Simulation with I = 200 and censoring rate of 50%

	Fam-simulation-300/Miss-20: Simulation with I = 300 and censoring rate of 20%
	Fam-simulation-300/Miss-50: Simulation with I = 300 and censoring rate of 50%

	Fam-simulation-500/Miss-20: Simulation with I = 500 and censoring rate of 20%
	Fam-simulation-500/Miss-50: Simulation with I = 500 and censoring rate of 50%

In Fam-simulation-300-sensitivity, we conduct simulations with I = 300 and censoring rate of 50%, but vary the prior hyper-parameters. 

In each folder, the main files are sim_1.R and sim_2.R, each of which runs 25 simulations for a total of 50 simulations. 

In the folder 'simulation analysis', the MCMC posterior samples in each simulation setting are stored in the corresponding folder (e.g., Fam-simulation-200-Miss-20). The main files are parameter_analysis.R, which reads in the posterior samples and produces the results for Table C.1 in Supplementary Materials, and penetrance_analysis.R and second_penetrance_analysis.R, which similarly produces the results for Table C.2 in Supplementary Materials. The numerical summaries are stored in the folder 'summary'. The file trace_plot_bw.R produces the trace plot as shown in Figure C.2 in Supplementary Materials.  

