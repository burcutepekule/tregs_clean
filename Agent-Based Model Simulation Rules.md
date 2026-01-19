Agent-Based Model Simulation Rules

## The Biological Story

Imagine the gut as a dynamic battleground where an intricate immune ballet unfolds. The epithelial layer, a single-cell-thick barrier separating our inner tissues from the microbial universe in our gut lumen, becomes injured in a localized area. This breach immediately triggers a cascade of danger signals called DAMPs (damage-associated molecular patterns) that radiate outward like distress beacons. As the barrier weakens, bacteria begin slipping through: some are harmless commensals that normally reside peacefully in the gut, while others are pathogenic invaders releasing their own toxic signals called PAMPs. Macrophages, the immune system's versatile sentinels, chemotax toward these combined danger signals and face a critical decision: should they become aggressive M1 warriors that produce tissue-damaging reactive oxygen species (ROS) to kill invaders, or should they adopt a gentler M2 resolving phenotype? This decision depends on the balance of danger signals (DAMPs plus PAMPs) versus safety signals (SAMPs). Meanwhile, regulatory T cells (Tregs) patrol the tissue, homing toward sites of damage. When they encounter macrophages presenting bacterial antigens, they must discriminate whether the antigen came from a harmless commensal or a dangerous pathogen. If Tregs perceive a commensal antigen, they activate and release SAMPs, anti-inflammatory signals that counter-balance the danger signals and help shift the local environment toward tolerance and repair. The macrophages are plastic, capable of switching between M0 resting, M1 inflammatory, and M2 resolving states based on which signals dominate their local environment. They maintain a "memory" or registry of recently engulfed bacteria that influences which antigens they present to Tregs. Critically, once a macrophage commits to a phenotype, it normally stays committed for a minimum period (active_age_limit) before reassessing its environment, providing phenotypic stability. However, the overwrite mechanism introduces flexibility: when enabled, if a strongly opposing signal appears, even young committed cells can switch phenotypes immediately, allowing rapid adaptation to sudden environmental changes. This represents the biological reality that while cells exhibit phenotypic commitment, sufficiently strong opposing signals can override this commitment before the typical reassessment period. The simulation captures this delicate balance between mounting an effective antimicrobial response to clear pathogens while avoiding excessive inflammation that causes collateral damage to the epithelium, ultimately determining whether the tissue successfully resolves the injury or descends into chronic inflammation.

1. General Parameters & Environment
Grid Size: 25×25 grid where each agent occupies one tile.
Time Steps: 5000 time steps.
Injury Site: An initial injury is defined as the middle 60% of the epithelial layer, representing a localized breach in the gut barrier.
Initial Populations: The simulation begins with 125 phagocytes (20% of grid surface), 125 Tregs (20% of grid surface), and 20 commensals in the lamina propria.
Signal Constraints: All signaling molecules (DAMPs, PAMPs, SAMPs, ROS) have maximum values of 1 to prevent unbounded accumulation. The combined danger signal (DAMPs + PAMPs) can reach a maximum of 2.
Stochasticity: Since the same parameter set can yield different outcomes depending on the random seed, each set is simulated 10 times. 
2. Epithelial Cells
Location: Epithelial cells are located at y=0 (effectively y=1 in the 1-indexed grid for interaction purposes).
Injury Levels: Epithelial cells have continuous injury levels starting from 0 (healthy), with a maximum constraint of max_level_injury.
Initial Injury: The simulation begins with epithelial cells within the injury site at injury level 1.
Injury Progression:
Pathogen-induced: Pathogens touching the basolateral side of epithelium (y=1) directly increase injury levels. Injury increment is calculated as: add_inj_pathogen = x0_in × (1 - exp(-k_in × pathogen_count)), using an exponential saturation function with maximum increment of x0_in. Commensals do not directly injure epithelial cells, but their detection at the basolateral side of the epithelium contribute to DAMPs release.
ROS-induced: ROS-induced injury is proportional to the amount by which local ROS exceeds the injury threshold: add_inj_ros = max(0, x0_in × (ros_mean - th_ROS_epith_injury)), with maximum increment of x0_in.
Maximum Injury: Injury level cannot exceed max_level_injury (typically 5).
Recovery: Injured epithelial cells recover stochastically with a certain probability (epith_recovery_chance) at each time step to reduce their injury level by 1.
DAMP Release:
Epithelial cells release DAMPs proportional to their level_injury.
Basolateral PRR stimulation by microbes (both commensals and pathogens) also triggers DAMP release from the epithelium using a logistic scaling function applied to the total microbe count at each epithelial position. The epithelium cannot differentiate between commensals and pathogens based on this initial PRR stimulation.
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
M2 (Anti-inflammatory/Resolving): Activated by safety signals. No ROS production, and identical engulfment activity to M1 (engulfment rate can be changed and sampled randomly, however, to be able to compare the scenarios more intuitively, we set the engulfment rate of M1 and M2 to the same value.)
Digestion:
Engulfment: At each time step, phagocytes attempt to engulf collocated microbes with probability = activity_engulf. Success rate depends on phenotype (M1 and M2 equal to each other for interpretability, but higher than M0).
Registry: Each phagocyte has a registry (interpreted as history or memory) of recently engulfed bacteria of a certain length of carrying capacity (cc_phagocyte). Entries are +1 for commensal, -1 for pathogen, 0 for empty. Registry is used for calculating the probabilities of presenting a commensal or pathogenic antigen on the MHC.
Digestion Time: Every digestion_time steps (default 1), the oldest entry in the registry is "digested" (removed), and new entries are added for newly engulfed microbes.
Activation and Deactivation Using Signal Dominance:
Danger Signal: Computed as danger_signal = DAMPs + PAMPs, representing the combined environmental threat.
Signal Difference Calculation: Phenotype decisions are based on normalized signal differences:
SAMPS_diff = max(0, avg_SAMPs - activation_threshold_SAMPs) / activation_threshold_SAMPs
DANGER_diff = max(0, danger_signal - activation_threshold_danger) / activation_threshold_danger
M0 → M1 Activation: If DANGER_diff > SAMPS_diff and DANGER_diff > 0, M0 phagocytes activate to M1.
M0 → M2 Prevention: Direct M0 → M2 transitions are prevented. This means Tregs can only suppress already-activated M1 macrophages (M1 → M2), not activate naive M0 macrophages. This biologically represents that tolerogenic signals require prior activation.
M1/M2 Reassessment: After reaching active_age_limit steps in their current phenotype, active phagocytes (M1 or M2) reassess their environment:
M1 → M2: If SAMPS_diff > DANGER_diff, M1 switches to M2 (SAMPs dominant).
M1 → M0: If both SAMPS_diff == 0 and DANGER_diff == 0, M1 deactivates to M0 (both signals below threshold).
M2 → M1: If DANGER_diff > SAMPS_diff, M2 switches to M1 (danger dominant).
M2 → M0: If both SAMPS_diff == 0 and DANGER_diff == 0, M2 deactivates to M0 (both signals below threshold).
Active Age Tracking: Each phagocyte tracks its active_age, incrementing each time step while in M1 or M2 phenotype. When phenotype switches, active_age resets to 1. M0 cells have active_age = 0.
Phenotype Plasticity with Overwrite Option:
Standard Mode (overwrite = 0): Phagocytes must remain in their committed phenotype (M1 or M2) for at least active_age_limit time steps before reassessing. This provides phenotypic stability and prevents rapid oscillations.
Overwrite Mode (overwrite = 1): Provides phenotypic flexibility by allowing "young" activated phagocytes (those with active_age < active_age_limit) to switch phenotypes immediately if the opposing signal becomes dominant, without waiting for active_age_limit. Specifically:
M1 → M2 early switch: If SAMPS_diff > DANGER_diff, even young M1 cells can switch to M2.
M2 → M1 early switch: If DANGER_diff > SAMPS_diff, even young M2 cells can switch to M1.
Biological Meaning: This represents the biological phenomenon where sufficiently strong opposing signals can overcome cellular commitment before the typical reassessment period. In the gut, this allows rapid response to sudden changes in the microbial or damage environment, such as a sudden influx of pathogens in a previously tolerogenic environment, or rapid resolution when Treg-derived SAMPs flood a previously inflammatory site.
ROS Production: Only M1 phagocytes produce ROS at their location. ROS diffuses after production.
Recruitment from Tissue Borders: Phagocytes can be recruited from bottom, left, and right borders (not epithelium) in response to local danger signals (DAMPs + PAMPs). Recruitment follows Poisson distribution with λ = recruitment_rate_danger × local danger signal. Newly recruited phagocytes start as M0 (resting) phenotype.
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
