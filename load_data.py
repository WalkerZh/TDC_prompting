import argparse
import json

import pandas as pd
from instruction_prompt import instruction_task_hub
from tdc_prompt import task_hub

SINGLE_ADME_TASK = ["Caco2_Wang", "Lipophilicity_AstraZeneca", "Solubility_AqSolDB",
                    "HydrationFreeEnergy_FreeSolv", "PPBR_AZ", "VDss_Lombardo", "Half_Life_Obach", "Clearance_Hepatocyte_AZ",
                    # ----- Regression
                    "PAMPA_NCATS", "HIA_Hou", "Pgp_Broccatelli", "Bioavailability_Ma", "BBB_Martins",
                    "CYP2C19_Veith", "CYP2D6_Veith", "CYP3A4_Veith", "CYP1A2_Veith", "CYP2C9_Veith",
                    "CYP2C9_Substrate_CarbonMangels", "CYP2D6_Substrate_CarbonMangels", "CYP3A4_Substrate_CarbonMangels"
                    # ----- Classification
                    ]
SINGLE_TOX_TASK = ["LD50_Zhu",
                    # ----- Regression
                   "hERG", "hERG_Karim", "AMES", "DILI", "Skin Reaction", "Carcinogens_Lagunin", "ClinTox",
                   # ----- Classification
                   "herg_central", "Tox21"
                   # ----- Multi subtasks
                   ]
SINGLE_HTS_TASK = ["SARSCoV2_Vitro_Touret", "SARSCoV2_3CLPro_Diamond", "HIV", 
                   "orexin1_receptor_butkiewicz", "m1_muscarinic_receptor_agonists_butkiewicz", "m1_muscarinic_receptor_antagonists_butkiewicz",
                   "potassium_ion_channel_kir2", "kcnq2_potassium_channel_butkiewicz", "cav3_t-type_calcium_channels_butkiewicz",
                   "choline_transporter_butkiewicz", "serine_threonine_kinase_33_butkiewicz", "tyrosyl-dna_phosphodiesterase_butkiewicz",
                   # ----- Classification
                   ]
CLASSIFICATION_TASK = ["PAMPA_NCATS", "HIA_Hou", "Pgp_Broccatelli", "Bioavailability_Ma", "BBB_Martins",
                    "CYP2C19_Veith", "CYP2D6_Veith", "CYP3A4_Veith", "CYP1A2_Veith", "CYP2C9_Veith",
                    "CYP2C9_Substrate_CarbonMangels", "CYP2D6_Substrate_CarbonMangels", "CYP3A4_Substrate_CarbonMangels",

                    "hERG", "hERG_Karim", "AMES", "DILI", "Skin Reaction", "Carcinogens_Lagunin", "ClinTox",
                    "herg_central", "Tox21",

                    "SARSCoV2_Vitro_Touret", "SARSCoV2_3CLPro_Diamond", "HIV", 
                    "orexin1_receptor_butkiewicz", "m1_muscarinic_receptor_agonists_butkiewicz", "m1_muscarinic_receptor_antagonists_butkiewicz",
                    "potassium_ion_channel_kir2", "kcnq2_potassium_channel_butkiewicz", "cav3_t-type_calcium_channels_butkiewicz",
                    "choline_transporter_butkiewicz", "serine_threonine_kinase_33_butkiewicz", "tyrosyl-dna_phosphodiesterase_butkiewicz",
                    ]
CLASSIFICATION_ADME_TOX_TASK = ["PAMPA_NCATS", "HIA_Hou", "Pgp_Broccatelli", "Bioavailability_Ma", "BBB_Martins",
                    "CYP2C19_Veith", "CYP2D6_Veith", "CYP3A4_Veith", "CYP1A2_Veith", "CYP2C9_Veith",
                    "CYP2C9_Substrate_CarbonMangels", "CYP2D6_Substrate_CarbonMangels", "CYP3A4_Substrate_CarbonMangels",

                    "hERG", "hERG_Karim", "AMES", "DILI", "Skin Reaction", "Carcinogens_Lagunin", "ClinTox",
                    "herg_central", "Tox21",
                    ]
CLASSIFICATION_HTS_TASK = ["SARSCoV2_Vitro_Touret", "SARSCoV2_3CLPro_Diamond", "HIV", 
                    "orexin1_receptor_butkiewicz", "m1_muscarinic_receptor_agonists_butkiewicz", "m1_muscarinic_receptor_antagonists_butkiewicz",
                    "potassium_ion_channel_kir2", "kcnq2_potassium_channel_butkiewicz", "cav3_t-type_calcium_channels_butkiewicz",
                    "choline_transporter_butkiewicz", "serine_threonine_kinase_33_butkiewicz", "tyrosyl-dna_phosphodiesterase_butkiewicz",
                    ]
REGRESSION_TASK = ["Caco2_Wang", "Lipophilicity_AstraZeneca", "Solubility_AqSolDB",
                "HydrationFreeEnergy_FreeSolv", "PPBR_AZ", "VDss_Lombardo", "Half_Life_Obach", "Clearance_Hepatocyte_AZ",

                "LD50_Zhu", "herg_central",
                ]

MULTI_DDI_TASK = ["DrugBank", "TWOSIDES"]
GENERATE_RETROSYN_TASK = ["USPTO"] # "USPTO-50K"

MOLECULENET_TASK = ["BBB_Martins", "ClinTox", "Tox21", "HIV"] # Bace / Sider

SINGLE_TASK = SINGLE_ADME_TASK + SINGLE_TOX_TASK + SINGLE_HTS_TASK
ALL_TASK = SINGLE_ADME_TASK + SINGLE_TOX_TASK + SINGLE_HTS_TASK + MULTI_DDI_TASK + GENERATE_RETROSYN_TASK

def output_file(dataset, outputs):
    '''
    mkdir / output
    '''
    return

def get_all_output(dataset, split, instruction_format=False, subtask=None, label_index=None):
    # splits = ["train", "valid", "test"]
    # ret = []
    
    # for ss in splits:
        # part_data = split[ss]
        # outputs = task_hub(dataset, part_data, subtask=subtask, label_index=label_index)
        # # for output in outputs[:2]:
        # #     print(output)
        # ret.extend(outputs)
    all_data = pd.concat([split['train'], split['valid'], split['test']], ignore_index=True)
    assert len(all_data) == len(split['train']) + len(split['valid']) + len(split['test'])
    if instruction_format:
        outputs = instruction_task_hub(dataset, all_data, subtask=subtask, label_index=label_index)
    else:
        outputs = task_hub(dataset, all_data, subtask=subtask, label_index=label_index)
    return outputs

def get_outputs_of_dataset(dataset, split_method, instruction_format=False, reg=1):
    if dataset in SINGLE_ADME_TASK:
        from tdc.single_pred import ADME
        data = ADME(name=dataset)
        split = data.get_split(method=split_method)
        # all_data = data.get_data()
        outputs = get_all_output(dataset, split, instruction_format=instruction_format)

    elif dataset in SINGLE_TOX_TASK:
        if dataset == "Tox21":
            from tdc.utils import retrieve_label_name_list
            label_list = retrieve_label_name_list(dataset)
            outputs = []
            for l in label_list:
                from tdc.single_pred import Tox
                data = Tox(name=dataset, label_name=l)
                split = data.get_split(method=split_method)
                output = get_all_output(dataset, split, instruction_format=instruction_format, subtask=l)
                outputs.extend(output)
        elif dataset == "herg_central":
            from tdc.utils import retrieve_label_name_list
            label_list = retrieve_label_name_list(dataset)
            outputs = []
            if reg == 0: # classification task
                from tdc.single_pred import Tox
                data = Tox(name=dataset, label_name=label_list[2])
                split = data.get_split(method=split_method)
                output = get_all_output(dataset, split, instruction_format=instruction_format, subtask=label_list[2])
                outputs.extend(output)
            else: # reg = 1 => regression task
                for l in label_list[:-1]:
                    from tdc.single_pred import Tox
                    data = Tox(name=dataset, label_name=l)
                    split = data.get_split(method=split_method)
                    output = get_all_output(dataset, split, instruction_format=instruction_format, subtask=l)
                    outputs.extend(output)

        else:
            from tdc.single_pred import Tox
            data = Tox(name=dataset)
            split = data.get_split(method=split_method)
            outputs = get_all_output(dataset, split, instruction_format=instruction_format)
    
    elif dataset in SINGLE_HTS_TASK:
        from tdc.single_pred import HTS
        data = HTS(name=dataset)
        split = data.get_split(method=split_method)
        outputs = get_all_output(dataset, split, instruction_format=instruction_format)
    
    elif dataset in MULTI_DDI_TASK:
        from tdc.multi_pred import DDI
        data = DDI(name=dataset)
        split = data.get_split(method=split_method)
        from tdc.utils import get_label_map
        if dataset == "DrugBank":
            label_index = get_label_map(name="DrugBank", task="DDI")
        elif dataset == "TWOSIDES":
            label_index = get_label_map(name="TWOSIDES", task="DDI", name_column="Side Effect Name")
        else:
            print("Dataset not exist!")
        outputs = get_all_output(dataset, split, instruction_format=instruction_format, label_index=label_index)

    elif dataset in GENERATE_RETROSYN_TASK:
        from tdc.generation import RetroSyn
        data = RetroSyn(name=dataset)
        split = data.get_split(method=split_method)
        outputs = get_all_output(dataset, split, instruction_format=instruction_format)
    
    else:
        print("Dataset not exist!")

    # print(f"Dataset {dataset} total output: {len(outputs)}.")
        
    return outputs

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", type=str, default="PAMPA_NCATS")
    parser.add_argument("--split", type=str, default="random")
    parser.add_argument("--instruction-format", type=bool, default=False)
    args = parser.parse_args()

    dataset = args.dataset
    if dataset == "Skin_Reaction":
        dataset = "Skin Reaction"
    
    outputs = get_outputs_of_dataset(dataset, args.split, instruction_format=args.instruction_format)

    print(f"Total examples of dataset \"{dataset}\": {len(outputs)}")
    print(f"Some examples of dataset \"{dataset}\":")
    for output in outputs[:3]:
        print(output)
    print(json.loads(outputs[0])['text'])
