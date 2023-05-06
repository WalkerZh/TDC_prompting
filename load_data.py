import argparse

import tdc_prompt


SINGLE_ADME_TASK = ["Caco2_Wang", "Lipophilicity_AstraZeneca", "Solubility_AqSolDB",
                    "HydrationFreeEnergy_FreeSolv", "PPBR_AZ", "VDss_Lombardo", "Half_Life_Obach", "Clearance_Hepatocyte_AZ",
                    "PAMPA_NCATS", "HIA_Hou", "Pgp_Broccatelli", "Bioavailability_Ma", "BBB_Martins",
                    "CYP2C19_Veith", "CYP2D6_Veith", "CYP3A4_Veith", "CYP1A2_Veith", "CYP2C9_Veith",
                    "CYP2C9_Substrate_CarbonMangels", "CYP2D6_Substrate_CarbonMangels", "CYP3A4_Substrate_CarbonMangels"
                    ]
SINGLE_TOX_TASK = ["LD50_Zhu",
                   "hERG", "hERG_Karim", "AMES", "DILI", "Skin Reaction", "Carcinogens_Lagunin", "ClinTox"
                   ]
SINGLE_HTS_TASK = []

def output_file(dataset, outputs):
    '''
    mkdir / output / (no-label?)
    '''
    return


parser = argparse.ArgumentParser()
parser.add_argument('--dataset', type=str, default="Caco-2")
parser.add_argument('--split', type=str, default="random")
args = parser.parse_args()

splits = ['train', 'valid', 'test']

if args.dataset == "Caco-2":
    from tdc.single_pred import ADME
    data = ADME(name='Caco2_Wang')
    split = data.get_split(method=args.split)
    for ss in splits:
        drug, y = split[ss]['Drug'], split[ss]['Y']
        outputs = tdc_prompt.get_Caco_2_prompt(drug, y)
