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

~~Total sample of 3 class: single(3,174,794), multi(255,000), generate(1,107,898)~~

#### Dataset_v3

###### Single-instance

classification: 2,410,684 can be obtained by `python load_all_data.py --task-class classification`

* ADME & TOX: 398,700 `python load_all_data.py --task-class classification_ADME_TOX`
* HTS: 2,011,984 ``python load_all_data.py --task-class classification_HTS``

regression: 641,529(613,786 from herg_central) (`<<reg_num_wei>>6.54321<<reg_num_wei>>`) can be obtained by `python load_all_data.py --task-class regression`

###### Multi-instance 

DDI: 255,110 (191637+63473) can be obtained by `python load_all_data.py --task-class multi`

###### Generation

retrosyn: 1,107,898 can be obtained by `python load_all_data.py --task-class generate`

## Single-instance Prediction Problem

#### ADME

- Caco2_Wang √ regression
- PAMPA_NCATS √ 2034 classification
- HIA_Hou √ 578 classification
- Pgp_Broccatelli √ 1218 classification
- Bioavailability_Ma √ 640 classification
- Lipophilicity_AstraZeneca √ regression
- Solubility_AqSolDB √ regression
- HydrationFreeEnergy_FreeSolv √ regression
- BBB_Martins × MoleculeNet
- PPBR_AZ √ regression

  **Note** : [Started from 0.3.7] this dataset contains assay across five species. Across species, the PPBR would behave rather differently even for the same drug. Thus, in default, it only returns the homo sapiens subset. If you would like to retrieve other species, you can use the following code:

  ```python
  data.get_other_species('Rattus norvegicus')
  # select from 'Canis lupus familiaris', 'Cavia porcellus', 'Homo sapiens', 'Mus musculus', 'Rattus norvegicus', 'all'
  ```
- VDss_Lombardo √ regression
- CYP2C19_Veith √ 12665 classification
- CYP2D6_Veith √ 13130 classification
- CYP3A4_Veith √ 12328 classification
- CYP1A2_Veith √ 12579 classification
- CYP2C9_Veith √ 12092 classification
- CYP2C9_Substrate_CarbonMangels √ 669 classification
- CYP2D6_Substrate_CarbonMangels √ 667 classification
- CYP3A4_Substrate_CarbonMangels √ 670 classification
- Half_Life_Obach √ regression
- Clearance_Hepatocyte_AZ √ regression

#### Tox

- LD50_Zhu √ regression
- hERG √ 655 classification
- herg_central √ classification & regression
  - Note:  hERG_inhib is a binary classification. Given a drug SMILES string, predict whether it blocks (1) or not blocks (0). This is **equivalent to whether hERG_at_10uM < -50**, i.e. whether the compound has an IC50 of less than 10µM.
- hERG_Karim √ 13445 classification
- AMES √ 7278 classification
- DILI √ 475 classification
- Skin_Reaction √ 404 classification
- Carcinogens_Lagunin √ 280 classification
- Tox21 × MoleculeNet
- Toxcast × MoleculeNet
- Clintox × MoleculeNet

#### HTS

- SARSCoV2_Vitro_Touret √ 1484 classification
- SARSCoV2_3CLPro_Diamond √ 880 classification
- HIV × MoleculeNet
- orexin1_receptor_butkiewicz √ 61829 classification
- m1_muscarinic_receptor_agonists_butkiewicz √ 61752 classification
- m1_muscarinic_receptor_antagonists_butkiewicz √ 61752 classification
- potassium_ion_channel_kir2 √ 301431 classification
- kcnq2_potassium_channel_butkiewicz √ 302343 classification
- cav3_t-type_calcium_channels_butkiewicz √ 100863 classification
- choline_transporter_butkiewicz √ 302242 classification
- serine_threonine_kinase_33_butkiewicz √ 319727 classification
- tyrosyl-dna_phosphodiesterase_butkiewicz √ 341324 classification

#### QM

> Not included, maybe add in later.

#### Yields

> Not included, maybe add in later.

## Multi-instance Prediction Problem

#### DDI

- DrugBank √(191,637)
- TWOSIDES √(63,473 pairs)
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
