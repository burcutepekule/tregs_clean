# AWS Compute Credits Application
## Computational Modeling of Regulatory T Cell Function in Gut Immune Homeostasis

---

## 1. DETAILED PROJECT DESCRIPTION

### Problem Statement

The mammalian intestine faces a unique immunological challenge: it must simultaneously maintain tolerance to trillions of beneficial commensal bacteria while mounting rapid defensive responses against pathogenic invaders. Regulatory T cells (Tregs) are critical mediators of this balance, but their role is paradoxical—in some contexts they prevent harmful inflammation and promote tissue healing, while in others they can suppress necessary immune responses and allow pathogen persistence. Understanding when and how Tregs help versus harm host defense is crucial for developing therapeutic interventions for inflammatory bowel disease (IBD), infectious diseases, and autoimmune disorders.

Despite extensive experimental work, we lack a mechanistic understanding of how local tissue conditions, microbial composition, and immune cell interactions determine whether Tregs will be beneficial or detrimental during gut injury and infection. The complexity of these multi-scale, dynamic interactions makes them difficult to study experimentally and necessitates computational approaches.

### Project Summary

This project employs spatially-explicit agent-based modeling (ABM) and complementary ordinary differential equation (ODE) models to systematically explore the parameter space governing immune responses in the gut lamina propria following epithelial barrier injury. The models simulate interactions between:

- **Epithelial cells** (25 cells with 6 injury states each)
- **Immune cells** (macrophages/dendritic cells with M0/M1/M2 phenotypes, regulatory T cells)
- **Microbes** (commensal bacteria and pathogens)
- **Signaling molecules** (DAMPs, PAMPs, SAMPs, ROS)

The ABM operates on a 25×25 spatial grid over 5,000 time steps, tracking individual cell movements, phenotype switching, and stochastic interactions. The model includes sophisticated biological mechanisms:

1. **Macrophage phenotype plasticity**: M0→M1 (pro-inflammatory) or M0→M2 (anti-inflammatory) based on local danger vs. safety signals
2. **Treg antigen discrimination**: Stochastic two-step process modeling antigen presentation and TCR recognition
3. **Spatial chemotaxis**: Cells migrate toward chemical gradients (DAMPs, PAMPs)
4. **Dynamic recruitment**: Immune cells recruited from tissue borders proportional to local inflammation
5. **ROS-mediated trade-offs**: Reactive oxygen species kill microbes but cause collateral tissue damage

Our central research questions are:
- Under what parameter regimes do Tregs accelerate vs. delay tissue healing?
- How does microbial composition (pathogen:commensal ratio) influence Treg impact?
- Can macrophage antigen discrimination compensate for absent or dysfunctional Tregs?
- What are the minimal mechanistic requirements for reproducing experimentally observed immune dynamics?

### Computational Challenge

The project requires comprehensive parameter space exploration using Latin Hypercube Sampling (LHS):

**Current status:**
- **87,504 parameter sets** generated (`lhs_parameters_della.csv`, 40 MB)
- **~12 scenarios** per parameter set (factorial combinations of: sterile/infectious injury, Tregs on/off, macrophage discrimination modes, ROS on/off)
- **10 replicates** per scenario to account for stochasticity
- **Total simulations needed**: ~10.5 million ABM runs

**Computational requirements per simulation:**
- ABM: 10-60 seconds per simulation (optimized with C++ via Rcpp, achieving 20-100× speedup)
- ODE: 1-2 seconds per simulation (deterministic approximation)
- Memory: ~500 MB per parallel job
- Output: ~50-100 KB per simulation (longitudinal time series data)

**Total computational estimate:**
- **ABM runs**: ~10.5M simulations × 30 sec average = 91,000 CPU-hours
- **ODE runs**: Additional ~5M simulations × 1.5 sec = 2,100 CPU-hours
- **Total**: ~93,000 CPU-hours of compute time
- **Storage**: ~750 GB for simulation outputs, ~200 GB for processed analysis

**Why we need AWS:**
Our current HPC allocation at Princeton (DELLA cluster) is shared and heavily utilized, with:
- Long queue times (24-48 hours typical wait)
- Limited burst capacity for parameter sweeps
- Storage constraints on scratch filesystems
- Difficulty scaling to full parameter space exploration

AWS will enable:
- Massive parallelization via EC2 Spot Instances (target: 1,000-2,000 concurrent cores)
- Flexible compute scaling to complete parameter sweeps within weeks rather than months
- Cost-effective storage with S3 for long-term data retention
- Integration with analysis pipelines using RStudio Server or Jupyter notebooks on EC2

### Proposed AWS Services

#### Primary Services

1. **Amazon EC2 (Elastic Compute Cloud)**
   - **Instance types**:
     - `c7i.8xlarge` (32 vCPUs, compute-optimized) for ABM simulations
     - `c7i.4xlarge` (16 vCPUs) for ODE simulations
     - `r7i.2xlarge` (8 vCPUs, memory-optimized) for data analysis and visualization
   - **Pricing strategy**: Spot Instances for 70-80% cost savings (simulations are interruptible and stateless)
   - **Estimated peak usage**: 500-1,000 vCPUs during parameter sweeps
   - **Duration**: 6-month intensive compute phase, then ongoing analysis

2. **Amazon S3 (Simple Storage Service)**
   - **Raw simulation outputs**: ~750 GB (RDS files with longitudinal time series)
   - **Processed datasets**: ~200 GB (aggregated metrics, classification results)
   - **Code repository mirror**: ~100 MB
   - **Versioning enabled** for reproducibility
   - **Lifecycle policies**: Transition infrequently accessed data to S3 Glacier after 6 months

3. **AWS Batch**
   - **Job orchestration**: Manage 87,504 parameter sets × scenarios as job arrays
   - **Automatic retries**: Handle Spot Instance interruptions gracefully
   - **Priority queues**: Balance exploratory vs. production runs
   - **Integration** with EC2 Spot Fleet for cost optimization

#### Supporting Services

4. **Amazon CloudWatch**
   - Monitor job completion rates
   - Track resource utilization and costs
   - Alert on failures or anomalies

5. **Amazon ECR (Elastic Container Registry)**
   - Store Docker images with R, Rcpp, deSolve, and model code
   - Enable reproducible compute environments
   - Version control for software dependencies

6. **AWS Lambda** (optional)
   - Trigger data processing pipelines when simulations complete
   - Automate aggregation of results from S3
   - Generate summary statistics and alerts

7. **Amazon RDS or DynamoDB**
   - Store simulation metadata (parameter set IDs, completion status, run times)
   - Enable efficient querying of completed vs. pending jobs
   - Track provenance for reproducibility

#### Analysis and Development

8. **Amazon EC2 (Interactive Instances)**
   - `r7i.2xlarge` for RStudio Server
   - Interactive data analysis, visualization, machine learning
   - Jupyter notebooks for documentation and tutorials

### Project Timeline and Key Milestones

**Month 1-2: Infrastructure Setup and Validation**
- Set up AWS Batch, S3 buckets, and ECR repositories
- Containerize R environment with all dependencies (Rcpp, deSolve, tidyverse)
- Implement job submission pipeline and monitoring dashboard
- Run pilot sweeps (1,000 parameter sets) to validate infrastructure
- Benchmark cost per simulation and optimize instance selection
- **Deliverable**: Validated compute pipeline, cost projection report

**Month 3-4: Primary Parameter Sweep (ABM)**
- Execute full ABM parameter sweep: 87,504 sets × 12 scenarios × 10 reps
- Continuous monitoring and quality control
- Intermediate data processing to identify interesting parameter regimes
- **Deliverable**: Complete ABM dataset (~750 GB), preliminary classification of parameter sets

**Month 5-6: ODE Model Calibration and Validation**
- Run ODE model on all parameter sets
- Calibration: Compare ABM vs. ODE predictions, tune ODE-specific parameters
- Identify parameter regimes where ODE approximation is valid vs. insufficient
- **Deliverable**: Validated ODE model, comparative analysis report

**Month 7-9: Data Analysis and Pattern Discovery**
- Machine learning classification: When do Tregs help/harm/have no effect?
- Parameter sensitivity analysis: Which biological parameters drive outcomes?
- Bifurcation analysis: Critical thresholds for Treg benefit vs. cost
- Comparison of macrophage discrimination vs. Treg-mediated tolerance
- **Deliverable**: Manuscript draft, interactive visualization dashboard

**Month 10-12: Targeted Exploration and Hypothesis Testing**
- Focused parameter sweeps around identified critical regimes
- Test specific hypotheses emerging from initial analysis
- Generate predictions for experimental validation
- Prepare educational materials and tutorials
- **Deliverable**: Submitted manuscript, public data release, tutorial documentation

---

## 2. PLAN FOR SHARING OUTCOMES

We are committed to full open science and will share all tools, data, and resources created during this project:

### Code and Software (Immediate Release)

1. **GitHub Repository** (already initiated, will be made public upon acceptance)
   - All simulation code (ABM and ODE implementations)
   - C++ acceleration functions (Rcpp)
   - Data generation and analysis scripts
   - Visualization tools and plotting functions
   - AWS deployment scripts (Dockerfiles, Batch job definitions, infrastructure-as-code)
   - **License**: MIT License for maximum reusability

2. **Docker Images** (via Docker Hub and Amazon ECR Public)
   - Production-ready containers with all dependencies
   - Tagged versions for reproducibility
   - Documentation for running simulations locally or on cloud platforms

3. **Documentation and Tutorials**
   - Comprehensive model documentation (already drafted: `Agent-Based Model Simulation Rules.md`, `ODE_README.md`)
   - Step-by-step tutorials for running simulations
   - Parameter interpretation guides
   - Case studies demonstrating how to test specific biological hypotheses

### Data (Public Release upon Publication)

4. **Simulation Outputs** (via Zenodo or Dryad)
   - Complete longitudinal datasets for all parameter sets and scenarios
   - Aggregated summary statistics and classification results
   - Metadata files with parameter mappings
   - **Size**: ~1 TB total (compressed)
   - **Format**: R-compatible RDS files + CSV exports for non-R users
   - **DOI-based citation** for academic credit

5. **Interactive Visualization Dashboard**
   - Web-based tool (via Shiny or Plotly Dash) for exploring parameter space
   - Hosted on AWS (if feasible) or Shiny server
   - Allow researchers to query outcomes without running simulations
   - Filter by parameter ranges, scenarios, outcomes

6. **Processed Datasets for Machine Learning**
   - Feature matrices for classification/regression tasks
   - Train/test splits for benchmarking
   - Enable method comparisons by other groups

### Educational Resources

7. **Preprint and Publication**
   - Preprint on bioRxiv immediately upon completion
   - Submission to peer-reviewed journal (e.g., PLoS Computational Biology, eLife, Cell Systems)
   - All data and code linked in manuscript

8. **Workshops and Presentations**
   - Tutorial workshops at conferences (e.g., Society for Mathematical Biology, American Society for Microbiology)
   - Webinars for immunology and computational biology communities
   - Lecture materials and recorded presentations

9. **Educational Modules**
   - Simplified versions of model for teaching computational immunology
   - Jupyter notebook tutorials suitable for graduate courses
   - Integration with existing educational platforms (e.g., nanoHUB, CoCalc)

### Principles
- **FAIR data**: Findable, Accessible, Interoperable, Reusable
- **Reproducibility**: All analyses fully documented and scripted
- **Community engagement**: Solicit feedback and contributions via GitHub
- **Accessibility**: Support for both R-fluent and non-R researchers

---

## 3. POTENTIAL FUTURE USE OF AWS BEYOND CREDIT AWARD DURATION

### By Individual Researcher

**Short-term (1-2 years post-award):**
- Refinement of model based on reviewer feedback and experimental validation
- Additional parameter sweeps for follow-up studies
- Integration of new biological mechanisms (e.g., T effector cells, IgA antibodies, epithelial regeneration dynamics)
- Sensitivity analysis and uncertainty quantification
- Personal cost likely sustainable via lab research funds for targeted analyses

**Long-term (3+ years):**
- Transition to hybrid model: local development + AWS for large-scale sweeps
- Use AWS for student training on cloud computing in computational biology
- Potential NIH R01 funding would support continued AWS usage for expanded model versions

### By Research Group (Metcalf Lab, Princeton EEB)

This project establishes infrastructure for the lab's broader computational immunology program:

1. **Related projects** that could leverage same AWS pipeline:
   - Within-host pathogen evolution under immune pressure
   - Age-structured immunity models requiring large parameter sweeps
   - Multi-scale models linking within-host dynamics to epidemiology

2. **Shared resources**:
   - Docker containers and Batch configurations reusable for other stochastic simulation projects
   - S3 data management practices applicable to empirical datasets (microbiome sequencing, flow cytometry)
   - CloudWatch monitoring templates for long-running computations

3. **Estimated ongoing usage**:
   - 1-2 intensive compute projects per year (each ~10,000-50,000 CPU-hours)
   - Lab budget allocation of $3,000-5,000/year for AWS (via grants)
   - Mix of on-demand for development and Spot for production runs

### By Broader Community

**Computational Immunology Community:**
- Model serves as template for other tissue-specific immune simulations (lung, skin, lymph node)
- Parameter exploration framework reusable for other agent-based models
- Benchmarking dataset for comparing ABM vs. ODE vs. PDE approaches

**Experimental Immunologists:**
- Interactive dashboard enables hypothesis generation without coding
- Predicted outcomes guide experimental design (e.g., which parameter regimes to test with mouse models)
- Potential for collaborative projects where experiments inform model refinement

**Education:**
- Medical and graduate schools teaching computational immunology
- Workshops at conferences generating episodic AWS usage
- Integration into existing online courses (e.g., Coursera, edX)

**Open Science Ecosystem:**
- Demonstrates best practices for reproducible, cloud-based computational research
- Contributes to NIH STRIDES initiative (cloud computing for biomedical research)
- Potential case study for AWS promotional materials on life sciences applications

### Sustainability Plan

1. **Grant funding**: Submit NIH R01 or NSF proposal highlighting AWS-based computational platform (success would provide 4-5 years of support)
2. **Cost optimization**: Transition mature, well-validated analyses to Spot Instances and reserved capacity for further savings
3. **Community contributions**: If tool gains traction, potential for user fees or institutional subscriptions to hosted dashboard
4. **Educational licensing**: Explore AWS Educate or research grants for continued student access

**Expected AWS cost trajectory:**
- **Year 1 (credit award)**: ~$30,000-40,000 equivalent in compute/storage
- **Year 2**: ~$8,000-12,000 (targeted analyses, grant-funded)
- **Year 3+**: ~$5,000/year baseline (maintenance, student projects, community support)

---

## 4. AWS EMPLOYEES CONTACTED

None at this time. We have not been in contact with any AWS employees prior to this application.

---

## 5. AWS PUBLIC DATASETS TO BE USED

While this project primarily generates new simulation data, we may integrate the following AWS Public Datasets for validation and contextualization:

### Primary Dataset

**NIH Human Microbiome Project (HMP)** (if available on AWS Open Data)
- **Purpose**: Validate commensal:pathogen ratios used in simulations
- **Usage**: Analyze gut microbiome composition distributions to inform realistic parameter ranges for microbial leakage rates and initial conditions
- **Integration**: Statistical distributions of bacterial abundances → model parameter priors

### Potential Secondary Datasets

**1000 Genomes Project** (if relevant genomic data becomes part of analysis)
- **Purpose**: If we extend model to include host genetic variation in immune parameters
- **Usage**: Population distributions of immune gene variants → parameter heterogeneity

**NIH Sequence Read Archive (SRA)** (if accessible via AWS)
- **Purpose**: Gut microbiome time series from IBD patients
- **Usage**: Compare model predictions of microbial dynamics during flares vs. remission to real patient data
- **Validation**: Test whether model can recapitulate observed pathogen bloom kinetics

**Cancer Genome Atlas (TCGA)** (indirectly relevant)
- **Purpose**: While focused on cancer, TCGA includes immune infiltrate data
- **Usage**: Inform macrophage:T cell ratios and immune cell densities in inflamed tissues

### Data Integration Strategy

1. **Parameter calibration**: Use empirical distributions from public data to constrain parameter ranges in LHS sampling
2. **Model validation**: Compare simulated outcomes to patient-derived time series when available
3. **Contextualization**: Demonstrate that model operates in biologically realistic regimes
4. **Publication**: Cite public datasets as validation sources, strengthening translational relevance

**Note**: If specific datasets listed above are not currently available on AWS Open Data, we will access them through standard repositories (NCBI, EBI) and upload to our S3 buckets for integration into analysis pipelines.

---

## SUMMARY

This project represents a significant computational immunology effort to understand the context-dependent roles of regulatory T cells in gut immune homeostasis. By leveraging AWS compute and storage resources, we will:

1. **Execute large-scale parameter sweeps** (93,000 CPU-hours) that are infeasible on shared HPC resources
2. **Develop open-source tools** for spatially-explicit immune simulation accessible to the broader community
3. **Generate publishable datasets** advancing our mechanistic understanding of immune tolerance vs. defense
4. **Establish cloud-based infrastructure** for ongoing computational immunology research
5. **Train students** in modern cloud computing for life sciences

All outputs will be made publicly available following FAIR principles, maximizing scientific impact and return on AWS's investment. The project has clear milestones, deliverables, and a sustainability plan for continued AWS usage beyond the credit award period.

We estimate requiring **$35,000-40,000 in AWS credits** over 12 months to complete this work.

---

## CONTACT INFORMATION

**Principal Investigator**: [Your Name]
**Institution**: [Your Institution]
**Department**: [Your Department]
**Email**: [Your Email]
**Project Website**: [To be created upon award]
**GitHub Repository**: https://github.com/[your-username]/tregs_clean (to be made public)

---

## SUPPORTING MATERIALS

Available upon request:
- Detailed technical specifications for AWS Batch job configurations
- Pilot data from preliminary parameter sweeps
- Letters of support from collaborating experimental immunologists
- Draft manuscript outline
- Institutional approval for cloud computing usage
- Data management plan detailing security and compliance measures

---

**Last Updated**: January 14, 2026
