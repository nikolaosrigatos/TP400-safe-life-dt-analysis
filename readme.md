The goal of this project is to perform a safe-life and damage-tolerance analysis of the TP400 propeller shaft, a component of the A400M military aircraft. Using the actual mean and alternating stresses recorded during various A400M mission profiles, three different approaches to safe-life and damage tolerance are carried out.

Task 1 focuses on a deterministic safe-life assessment using the “top of the scatter” approach. The material S–N curve is used to fit a log-likelihood function, while the conservative material curve (mu minus three sigma) and the extreme equivalent stress curve (mu plus three sigma) are used to calculate fatigue damage through Miner’s rule. The component life is then estimated as the inverse of the accumulated damage. A probabilistic version of Miner’s rule is also applied using a Monte Carlo analysis, demonstrating that the mean probabilistic result matches the deterministic top-of-the-scatter prediction.

Task 2 begins with an elastoplastic fracture assessment to determine the critical crack size. Two crack configurations are considered: a circumferential through-wall crack and a plate with a central crack. After assuming an initial maximum crack size, the maximum stress intensity factor and maximum ligament-yielding parameter are calculated for each configuration. A correction factor for ligament yielding is applied to determine the stress intensity factor for different crack sizes. The critical crack size for each configuration is the smallest value at which either fracture instability or plastic collapse occurs. The governing critical crack size for the component is the smallest among the two configurations.

A deterministic crack-growth analysis is then performed using the mission stress profiles, with the loading sequence randomized. Crack growth (a vs. N) is calculated using NASGRO. Based on the most conservative a-vs-N curve, an inspection plan is developed to achieve a 3000-hour inspection interval. Three non-destructive inspection methods are considered: eddy current, liquid penetrant with directed inspection, and liquid penetrant with full inspection. For each method, the inspection interval and required number of inspections are determined.

Task 3 involves a probabilistic damage-tolerance assessment. Starting from the anomalies distribution for circular holes given in AC 33-70.2, a probabilistic analysis is performed considering the presence of eight holes on each of the four shafts and applying a weakest-link approach. Finally, a Monte Carlo–based inspection plan is generated following the scheme proposed in AC 33-70.2 (Figure A7-6).

- Code contains all the MATLAB files
  - FAD_config1.m calculates the Failure Assessment Diagram (FAD) for configuration one (shown in TP400_Nikolaos_Rigatos.pdf, slide 18)
  - FAD_config2.m calculates the Failure Assessment Diagram (FAD) for configuration two (shown in TP400_Nikolaos_Rigatos.pdf, slide 18)
  - MC_inspection_plan.m performs a Monte Carlo (MC) analysis to produce a probabilistic inspection plan
  - PND_MC.m and PND_calc.m calculate the Probability of Non-Detection (PND) for the MC analysis and the deterministic damage-tolerance analysis
  - TP400.m is the main file that runs all functions
  - damage_calc.m performs a top-of-the-scatter life calculation
  - inspection_plan_profile1.m calculates the inspections required to achieve a PND of less than 2e-5
  - interpolate_curves.m interpolates the crack-versus-cycles graphs extracted from NASGRO
  - loglik_3par_norm_runouts.m calculates the log-likelihood function for the S-N curve
  - prob_assessment.m performs a probabilistic top-of-the-scatter assessment
  - prob_damage_tol.m performs the probabilistic damage-tolerance analysis for a probability of exceedance of 2e-5
  - write_SF_files.m writes the NASGRO files required for the deterministic damage-tolerance analysis
    
- Fatigue Data contains the Excel files for the S-N curve
- Mission Profiles contains the mean and alternating stress data for each mission profile
- Monte Carlo contains the crack-versus-cycles data for different initial crack lengths, used for interpolation via interpolate_curves
- NAGRO contains the NASGRO results for Task 2
- PDF contains all PDFs used for bibliography, tasks, and references
- POD Curves contains the probability-of-detection data for each inspection tool
- Results contains the .png files of all results
- SF_NASGRO contains the NASGRO input files for Task 2
- TP_400_Nikolaos_Rigatos.pdf and TP_400_Nikolaos_Rigatos.pptx are the short presentations of the project
