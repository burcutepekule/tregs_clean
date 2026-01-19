Agent-Based Model Simulation Rules

1. General Parameters & Environment
Grid Size: 25×25 grid where each agent occupies one tile.
Time Steps: 5000 time steps.
Injury Site: An initial injury is defined as the middle 60% of the epithelial layer, representing a localized breach in the gut barrier.
Initial Populations: The simulation begins with 125 phagocytes (20% of grid surface), 125 Tregs (20% of grid surface), and 20 commensals in the lamina propria.
Signal Constraints: All signaling molecules (DAMPs, PAMPs, SAMPs, ROS) have maximum values of 1 to prevent unbounded accumulation.
Stochasticity: Since the same parameter set can yield different outcomes depending on the random seed, each set is simulated 10 times. 
2. Epithelial Cells
Location: Epithelial cells are located at y=0 (effectively y=1 in the 1-indexed grid for interaction purposes).
Injury Levels: Epithelial cells have 6 states:
0: Healthy
1 to 5: Injured, with 5 being the maximum injury level.
Initial Injury: The simulation begins with epithelial cells within the injury site at injury level 1.
Injury Progression:
Pathogen-induced: Pathogens touching the basolateral side of epithelium (y=1) directly increase injury levels. Injury is calculated using a logistic function that maps pathogen count to injury level (0-5). Commensals do not directly injure epithelial cells, but their detection at the basolateral side of the epithelium contribute to DAMPs release.
ROS-induced: If the mean ROS level in the vicinity of an epithelial cell exceeds a certain threshold (th_ROS_epith_injury), its injury level increases by 1.
Maximum Injury: Injury level cannot exceed 5 (max_level_injury).
Recovery: Injured epithelial cells recover stochastically with a certain probability (epith_recovery_chance) at each time step to reduce their injury level by 1.
DAMP Release:
Epithelial cells release DAMPs proportional to their level_injury.
Basolateral PRR stimulation by microbes (both commensals and pathogens) also triggers DAMP release from the epithelium. The epithelium cannot differentiate between commensals and pathogens based on this initial PRR stimulation.
3. Microbes
Types: Commensals (non-pathogenic bacteria that normally colonize the gut) and Pathogens (virulent bacteria that damage epithelial cells and produce PAMPs).
Movement: Microbes move randomly. At the epithelium (y=1), they can only move deeper into tissue or laterally. Within the lamina propria (y>1), they can move in any of 8 directions or stay still.
Leakage into Lamina Propria:
Pathogens: Leak into the lamina propria from injured epithelial cells, with rate proportional to average injury level and injury site length. If it's a sterile injury, no pathogens leak.
Commensals: Leak at a baseline rate across the entire epithelium during homeostasis and no injured epithelial cells, and at an increased rate from injured epithelial cells.
ROS-induced Death: Microbes die if the local ROS concentration at their location exceeds a certain threshold (th_ROS_microbe).
4. Phagocytes (Macrophages/Dendritic Cells)
Population: Initial number of 125 phagocytes (20% of grid x grid). Phagocytes do not proliferate or die during simulation but can be dynamically recruited from tissue borders in response to danger signals.
Movement: Phagocytes chemotax toward combined danger signals (DAMPs + PAMPs). If no gradient exists, they move randomly.
Phenotypes: Phagocytes can have three phenotypes:
M0 (Resting): Default state with no ROS production and low engulfment activity.
M1 (Pro-inflammatory): Activated by danger signals. High ROS production and high engulfment activity. 
M2 (Anti-inflammatory/Resolving): Activated by safety signals. No ROS production, and identical engulfment activity to M2 (engulfment rate can be changed and sampled randomly, however, to be able to compare the scenarios more intuitively, we set the engulfment rate of M1 and M2 to the same value.)
Digestion:
Engulfment: At each time step, phagocytes attempt to engulf collocated microbes with probability = activity_engulf. Success rate depends on phenotype (M1 and M2 equal to each other for interpretability, but higher than M0).
Registry: Each phagocyte has a registry (interpreted as history or memory) of recently engulfed bacteria of a certain length of carrying capacity (cc_phagocyte). Entries are +1 for commensal, -1 for pathogen, 0 for empty. Registry is used for calculating the probabilities of presenting a commensal or pathogenic antigen on the MHC. 
Digestion Time: Every digestion_time steps (default 1), the oldest entry in the registry is "digested" (removed), and new entries are added for newly engulfed microbes. 
Activation and Deactivation: After active_age_limit steps, active phagocytes (M1 or M2) reassess their environment:
If avg_DAMPs >= activation_threshold_DAMPs and avg_DAMPs > avg_SAMPs, the phagocyte becomes (or remains) M1.
If avg_SAMPs >= activation_threshold_SAMPs and avg_SAMPs > avg_DAMPs, the phagocyte becomes (or remains) M2. This means an M1 phagocyte can differentiate into an M2 if SAMPs become the dominant signal.
If both avg_SAMPs < activation_threshold_SAMPs and avg_DAMPs < activation_threshold_DAMPs, the phagocyte returns to the M0 (resting) state.
Phenotype Plasticity (M1 ↔ M2 ↔ M0): When any cell (including phagocytes and T cells) switch their phenotype, they stay in that phenotype for active_age_limit time steps no matter the stimulus. After active_age_limit time steps, M1 and M2 phagocytes reassess their environment and can switch phenotypes based on signal dominance.
ROS Production: Only M1 phagocytes produce ROS at their location. ROS diffuses after production.
Recruitment from Tissue Borders: Phagocytes can be recruited from top, left, and right borders (not epithelium) in response to local danger signals (DAMPs + PAMPs). Recruitment follows Poisson distribution with λ = recruitment_rate_danger × local danger signal.
Macrophage Specificity (optional): When macspec_on > 0, phagocytes use both environmental signals AND engulfed bacteria for polarization decisions. Assuming macrophages also have receptors like Tregs, they first identify the antigen with accuracy proportional to mac_discrimination_efficiency. They then use a danger-biased logic: M1 activation uses OR logic (either environmental danger OR perceived pathogen is sufficient), while M2 activation uses AND logic (requires BOTH environmental safety AND perceived commensal engulfment).

5. Tregs (Regulatory T cells)
Population: A fixed number of 125 Tregs (20% of grid x grid). Tregs do not proliferate, die, or are recruited during simulation.
Movement: Tregs chemotax toward DAMPs (sites of tissue damage). If no DAMP gradient exists or randomize_tregs = 1, they move randomly.
Phenotypes: Tregs can be:
Resting (phenotype = 0): Default state, no SAMP production.
Activated (phenotype = 1): Produces SAMPs at their location.
Activation:
If allow_tregs_to_do_their_job is TRUE: Tregs activate if they are in the vicinity (treg_vicinity_effect = 1) of an M1 or M2 phagocyte that is presenting an antigen perceived as a commensal by the Treg (this can be achieved by either the phagocyte presenting a commensal antigen and that being correctly perceived by the Treg, or phagocyte presenting a pathogenic antigen and that being erroneously perceived as a commensal antigen by the Treg)

Stochastic Two-Step Activation Process:
Step 1 - Antigen presentation: Stochastic draw determines which antigen type (commensal or pathogen) is presented to the Treg, based on actual bacterial composition in the phagocyte's registry.
Step 2 - Antigen identification: Treg identifies the presented antigen with accuracy proportional to treg_discrimination_efficiency. High efficiency means accurate identification; low efficiency leads to misidentification.
Activation decision: Treg activates if it perceives a commensal antigen (regardless of actual truth). All nearby Tregs within treg_vicinity_effect activate simultaneously.
SAMP Production: Activated Tregs produce SAMPs at their location.
Deactivation: After active_age_limit time steps in activated state, Tregs return to resting and stop SAMP production unless activated again

6. Signaling Molecules (DAMPs, SAMPs, ROS)

DAMPs (Damage-Associated Molecular Patterns):
Sources: Injured epithelial cells (proportional to injury level) and microbes at epithelium triggering PRR stimulation.
Function: Signals tissue damage and danger, activates phagocytes → M1, attracts phagocytes and Tregs.
PAMPs (Pathogen-Associated Molecular Patterns):
Sources: Represents toxins such as LPS, flagellin, and peptidoglycans, etc. released by pathogens at their location, proportional to pathogen count. 
Function: Signals pathogen presence. When combined with DAMPs, it forms danger signals for phagocyte activation and attraction. 
SAMPs (Safety/Suppression-Associated Molecular Patterns):
Sources: Represents anti-inflammatory cytokines such as IL-10 and TGF-β released by activated Tregs.
Function: Signals immunosuppression and tissue repair, activates phagocytes → M2, counter-balances danger signals. It is important to note that even if they are not high enough to shift phagocytes to an M2 phenotype, they can still prevent them to become M1 by competing with the danger signals.
ROS (Reactive Oxygen Species): 
Sources: Released by M1 phagocytes.
Function: Antimicrobial activity (kills microbes) and causes epithelial injury (collateral damage).
Control scenario: When control = 1, ROS production is disabled to test if infection can resolve without ROS. If this is the case, the corresponding parameter combination is discarded since this parameter set is not observed in nature (otherwise why would you need ROS and all that damage anyway if you could resolve everything with an anti-inflammatory phenotype?)
Diffusion: All signals diffuse using a discrete Laplacian operator with reflective boundary conditions. Linear decay at each time step prevents unbounded accumulation. All signals capped at maximum value of 1, meaning that danger_signals = PAMPs+DAMPs is capped at 2. 
