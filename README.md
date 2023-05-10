# TDC_prompting

This is a repo for building TDC(Therapeutics Data Commons) data prompts.

[toc]

## Usage

#### Installation

```bash
# Install necessary packages
pip install rdkit==2022.9.5 
pip install PyTDC # https://github.com/mims-harvard/TDC

# clone the repo
git clone https://github.com/WalkerZh/TDC_prompting.git
cd TDC_prompting
```

#### Load data

```bash
mkdir -p prompt_data
python load_all_data.py	--output-dir prompt_data	# by doing this, you can get a complete prompting data in folder prompt_data

# Or you can use command below to get one specified prompting dataset
python load_data.py --dataset <dataset_name>
```

All prompts are saved in `tdc_prompt.py`

## Description

TDC covers a wide range of therapeutics tasks with varying data structures. Thus, we organize it into three layers of hierarchies. First, we divide into three distinctive machine learning **problems**:

- Single-instance prediction `single_pred`: Prediction of property given individual biomedical entity.
- Multi-instance prediction `multi_pred`: Prediction of property given multiple biomedical entities.
- Generation `generation`: Generation of new biomedical entity.

## Single-instance Prediction Problem

#### ADME

- Caco2_Wang √
- PAMPA_NCATS √
- HIA_Hou √
- Pgp_Broccatelli √
- Bioavailability_Ma √
- Lipophilicity_AstraZeneca √
- Solubility_AqSolDB √
- HydrationFreeEnergy_FreeSolv √
- BBB_Martins √
- PPBR_AZ √

  **Note** : [Started from 0.3.7] this dataset contains assay across five species. Across species, the PPBR would behave rather differently even for the same drug. Thus, in default, it only returns the homo sapiens subset. If you would like to retrieve other species, you can use the following code:

  ```python
  data.get_other_species('Rattus norvegicus')
  # select from 'Canis lupus familiaris', 'Cavia porcellus', 'Homo sapiens', 'Mus musculus', 'Rattus norvegicus', 'all'
  ```
- VDss_Lombardo √
- CYP2C19_Veith √
- CYP2D6_Veith √
- CYP3A4_Veith √
- CYP1A2_Veith √
- CYP2C9_Veith √
- CYP2C9_Substrate_CarbonMangels √
- CYP2D6_Substrate_CarbonMangels √
- CYP3A4_Substrate_CarbonMangels √
- Half_Life_Obach √
- Clearance_Hepatocyte_AZ √

#### Tox

- LD50_Zhu √
- hERG √
- herg_central √
- hERG_Karim √
- AMES √
- DILI √
- Skin_Reaction √
- Carcinogens_Lagunin √
- Tox21 √
- Toxcast ×
- Clintox √

#### HTS

- SARSCoV2_Vitro_Touret √
- SARSCoV2_3CLPro_Diamond √
- HIV √
- orexin1_receptor_butkiewicz √
- m1_muscarinic_receptor_agonists_butkiewicz √
- m1_muscarinic_receptor_antagonists_butkiewicz √
- potassium_ion_channel_kir2 √
- kcnq2_potassium_channel_butkiewicz √
- cav3_t-type_calcium_channels_butkiewicz √
- choline_transporter_butkiewicz √
- serine_threonine_kinase_33_butkiewicz √
- tyrosyl-dna_phosphodiesterase_butkiewicz √

#### QM

> Not included, maybe add in later.

#### Yields

> Not included, maybe add in later.

## Multi-instance Prediction Problem

#### DDI

- DrugBank √
- TWOSIDES √
  - **Note**: Some side-effects may be a verb (e.g. agitated), making the prompt not clear or not correct

## Generation Problem

#### RetroSyn

- USPTO-50K √

  **Note:** To get the reaction types of each reaction, you can type:

  ```python
  from tdc.utils import get_reaction_type
  get_reaction_type('USPTO-50K')
  ```
- USPTO √
