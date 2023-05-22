import json
import os
import re

from rdkit import Chem  # rdkit==2022.9.5
from tqdm import tqdm

# format ref
# {"text": "We can conclude that the screening result of ability to inhibit HIV replication of <<|mol0|>> is inactive .", "entities": {"<<|mol0|>>": {"smiles": "CCOP(=O)(Nc1cccc(Cl)c1)OCC"}}}

system_message = "Below is an instruction that describes a task. Write a response that appropriately completes the request."

DRUG_Y = {
    "Caco2_Wang": {
        "common": r"The experimental result on the rate of {smiles} passing through the Caco-2 cells is <<reg_num_wei>>{label}<<reg_num_wei>>.",
    }, # regression log cm/s
    "PAMPA_NCATS": {
        "common": r"The compound {smiles} has {label} permeability in parallel artificial membrane permeability assay(PAMPA assay).",
        "instruction_version": {
            "instruction": "How permeable is the drug represented by SMILES below to the cellular membrane in the parallel artificial membrane permeability assay(PAMPA assay)?", 
            # Does the drug have high permeability or low-to-moderate permeability across the cellular membrane?
            "input": r"{smiles}",
            "positive_response": "It has high permeability across the cellular membrane.",
            "negative_response": "It has low-to-moderate permeability across the cellular membrane.",
        }
    },
    "HIA_Hou": {
        "common": r"The HIA(Human Intestinal Absorption) activity of {smiles} is {label}.",
        "instruction_version": {
            "instruction": "Can the drug represented by SMILES below be absorbed from the human gastrointestinal system into the bloodstream of the human body?",
            "input": r"{smiles}",
            "positive_response": "Yes, it can be absorbed from...", # TODO
            "negative_response": "No, it can't be absorbed.",
        },
    },
    "Pgp_Broccatelli": {
        "common": r"The drug {smiles} is {label} in inhibiting P-glycoprotein(Pgp).",
        "instruction_version": {
            "instruction": "Can the drug represented by SMILES below inhibit P-glycoprotein (Pgp)?", # detailed description?
            "input": r"{smiles}", # smiles
            "positive_response": "Yes, it can inhibit P-glycoprotein.",
            "negative_response": "No, it can't inhibit P-glycoprotein.",
        },
    },
    "Bioavailability_Ma": {
        "common": r"The oral bioavailability activity of {smiles} is {label}.",
        "instruction_version": {
            "instruction": "How is the oral bioavailability of this drug represented by SMILES below?",
            "input": r"{smiles}", # smiles
            "positive_response": "The oral bioavailability of it is active.", # or It is active?
            "negative_response": "The oral bioavailability of it is inactive.",
        },
    },

    "Lipophilicity_AstraZeneca": r"The ability of {smiles} to dissolve in a lipid (e.g. fats, oils) environment is <<reg_num_wei>>{label}<<reg_num_wei>>.", # regression # units? # experimental results of octanol/water distribution coefficient (logD at pH 7.4)
    "Solubility_AqSolDB": r"The log solubility of {smiles} in water is <<reg_num_wei>>{label}<<reg_num_wei>> in mols per litre.", # regression # units: log solubility in mols per litre
    # "Solubility_AqSolDB": r"The ability of {smiles} to dissolve in water is <<reg_num_wei>>{label}<<reg_num_wei>>.", # regression # units? # measured log solubility in mols per litre
    "HydrationFreeEnergy_FreeSolv": r"The log hydration free energy of {smiles} in water is <<reg_num_wei>>{label}<<reg_num_wei>> in kcals per mol.", # regression # units? # kcal/mol
    
    "BBB_Martins": {
        "common": r"The ability of {smiles} to penetrate the Blood-Brain Barrier(BBB) is {label}.",
        "instruction_version": {
            "instruction": "Can the drug represented by SMILES below penetrate the blood-brain barrier(BBB)?",
            "input": r"", # smiles
            "positive_response": "Yes, it can penetrate the blood-brain barrier.",
            "negative_response": "No, it can't penetrate the blood-brain barrier.",
        },
    },

    "PPBR_AZ": r"The percentage of {smiles} bound to plasma proteins in the blood, i.e. its human plasma protein binding rate (PPBR) is <<reg_num_wei>>{label}<<reg_num_wei>>%.", # regression # units %
    "VDss_Lombardo": r"The volume of distribution at steady state (VDss) of {smiles}, i.e. the degree of its concentration in body tissue compared to concentration in blood is <<reg_num_wei>>{label}<<reg_num_wei>>.", # regression # units? # improve
    
    "CYP2C19_Veith": {
        "common": r"The drug {smiles} is {label} in inhibiting the CYP2C19 gene, which provides instructions for making an enzyme called the endoplasmic reticulum.",
        "instruction_version": {
            "instruction": "How is the activity of the drug represented by SMILES below in inhibiting the CYP2C19 gene?",
            # Is the drug represented by SMILES below active in inhibiting the CYP2C19 gene?
            "input": r"", # smiles
            "positive_response": "",
            "negative_response": "",
        },
    },

    "CYP2D6_Veith": r"The drug {smiles} is {label} in inhibiting the CYP2D6, which is primarily expressed in the liver.",
    "CYP3A4_Veith": r"The drug {smiles} is {label} in inhibiting the CYP3A4, which is an important enzyme that oxidizes small foreign organic molecules in the body.",
    "CYP1A2_Veith": r"The drug {smiles} is {label} in inhibiting the CYP1A2, which is able to metabolize some PAHs to carcinogenic intermediates.",
    "CYP2C9_Veith": r"The drug {smiles} is {label} in inhibiting the CYP2C9, which plays a major role in the oxidation of both xenobiotic and endogenous compounds.",
    "CYP2C9_Substrate_CarbonMangels": {
        "common": r"The drug {smiles} {label} a substrate to the CYP2C9 enzyme.",
        "instruction_version": {
            "instruction": "Is the drug represented by SMILES below a substrate to the CYP2C9 enzyme?",
            "input": r"", # smiles
            "positive_response": "Yes, it is.",
            "negative_response": "No, it is not.",
        },
    },
    "CYP2D6_Substrate_CarbonMangels": r"The drug {smiles} {label} a substrate to the CYP2D6 enzyme.",
    "CYP3A4_Substrate_CarbonMangels": r"The drug {smiles} {label} a substrate to the CYP3A4 enzyme.",

    "Half_Life_Obach": r"The half life duration of {smiles} is <<reg_num_wei>>{label}<<reg_num_wei>> hours.", # regression # units hours
    "Clearance_Hepatocyte_AZ": r"The volume of plasma cleared of {smiles} over a specified time period, i.e. its drug clearance is <<reg_num_wei>>{label}<<reg_num_wei>>.", # regression # units?
    # ------ ADME
    "LD50_Zhu": r"The acute toxicity of {smiles}, i.e. its most conservative dose that can lead to lethal adverse effects is <<reg_num_wei>>{label}<<reg_num_wei>>.", # regression # units?
    "hERG": {
        "common": r"The drug {smiles} {label} the Human ether-à-go-go related gene(hERG), which is crucial for the coordination of the heart's beating.", # the type of Y is 'float' although it's a classification task
        "instruction_version": {
            "instruction": "Does the drug represented by SMILES below block the Human ether-à-go-go related gene(hERG)?",
            "input": r"", # smiles
            "positive_response": "Yes, it blocks hERG.",
            "negative_response": "No, it does not block hERG.",
        },
    },
    "hERG_Karim": r"The drug {smiles} {label} the Human ether-à-go-go related gene(hERG), which is crucial for the coordination of the heart's beating.",
    "AMES": {
        "common": r"Mutagenicity means the ability of a drug to induce genetic alterations, the drug {smiles} {label} mutagenic.",
        "instruction_version": {
            "instruction": "Is the drug represented by SMILES below mutagenic?",
            "input": r"", # smiles
            "positive_response": "Yes, it is mutagenic.",
            "negative_response": "No, it's not mutagenic.",
        },
    },
    "DILI": {
        "common": r"The drug {smiles} {label} cause liver injury.",
        "instruction_version": {
            "instruction": "Please tell me whether the drug represented by SMILES below can cause liver injury or not.",
            "input": r"", # smiles
            "positive_response": "This drug can cause liver injury.",
            "negative_response": "This drug cannot cause liver injury.",
        },
    },
    "Skin Reaction": {
        "common": r"The drug {smiles} {label} cause skin reaction.",
        "instruction_version": {
            "instruction": "Can the drug represented by SMILES below cause skin reaction?",
            "input": r"", # smiles
            "positive_response": "Yes, This drug can cause skin reaction.",
            "negative_response": "No, This drug cannot cause skin reaction.",
        },
    },
    "Carcinogens_Lagunin": {
        "common": r"The drug {smiles} {label} cause carcinogen, which is any substance, radionuclide, or radiation that promotes carcinogenesis, the formation of cancer.",
        "instruction_version": {
            "instruction": "Can the drug represented by SMILES below cause carcinogen?",
            "input": r"", # smiles
            "positive_response": "Yes, This drug can cause carcinogen.",
            "negative_response": "No, This drug cannot cause carcinogen.",
        },
    },
    "ClinTox": r"The drug {smiles} has {label} clinical trials for toxicity reasons.",
    #  ------ Tox
    "SARSCoV2_Vitro_Touret": {
        "common": r"The activity of {smiles} against SARS-CoV-2 is {label}.",
        "instruction_version": {
            "instruction": "How is the activity of the drug represented by SMILES below against SARS-CoV-2?",
            # Is the drug represented by SMILES below active against SARS-CoV-2?
            "input": r"", # smiles
            "positive_response": "",
            "negative_response": "",
        },
    },
    "SARSCoV2_3CLPro_Diamond": r"The activity of {smiles} against SARS-CoV-2 3C-like protease(SARSCoV2 3CL) is {label}.",
    "HIV": r"The ability of {smiles} in inhibiting HIV replication, i.e. its activity against HIV virus is {label}.", # MoleculeNet
    "orexin1_receptor_butkiewicz": r"The activity of {smiles} against Orexin1 Receptor(HCRTR1) is {label}.",
    "m1_muscarinic_receptor_agonists_butkiewicz": r"The activity of {smiles} against M1 muscarinic receptor agonists is {label}.",
    "m1_muscarinic_receptor_antagonists_butkiewicz": r"The activity of {smiles} against M1 muscarinic receptor antagonists is {label}.",
    "potassium_ion_channel_kir2": r"The activity of {smiles} against potassium ion channel Kir2.1(KCNJ2) is {label}.",
    "kcnq2_potassium_channel_butkiewicz": r"The activity of {smiles} against KCNQ2 potassium channel(Kv7.2) is {label}.",
    "cav3_t-type_calcium_channels_butkiewicz": r"The activity of {smiles} against Cav3 T-type calcium channels is {label}.",
    "choline_transporter_butkiewicz": r"The activity of {smiles} against Choline transporter is {label}.", # , a protein that mediates the transport of choline across cell membranes,
    "serine_threonine_kinase_33_butkiewicz": r"The activity of {smiles} against Serine/threonine kinase 33(STK33) is {label}.",
    "tyrosyl-dna_phosphodiesterase_butkiewicz": r"The activity of {smiles} against Tyrosyl-DNA phosphodiesterase 1(TDP1) is {label}.",
    # ------ HTS
}

HERG_CENTRAL = { # DRUG_Y format
    "hERG_at_1uM": r"The percent inhibition of {smiles} to the Human ether-à-go-go related gene(hERG) at a 1µM concentration is <<reg_num_wei>>{label}<<reg_num_wei>>.", # regression # units?
    "hERG_at_10uM": r"The percent inhibition of {smiles} to the Human ether-à-go-go related gene(hERG) at a 10µM concentration is <<reg_num_wei>>{label}<<reg_num_wei>>.", # regression # units?
    "hERG_inhib": r"The drug {smiles} {label} the Human ether-à-go-go related gene(hERG), which is crucial for the coordination of the heart's beating.", # whether hERG_at_10uM < -50, i.e. whether the compound has an IC50 of less than 10µM.
}

TOX21 = { # DRUG_Y format
    "NR-AR": {
        "common": r"The compound {smiles} {label} activate the Androgen Receptor(AR).",
        "instruction_version": {
            "instruction": "Can the compound represented by SMILES below activate the Androgen Receptor(AR)?",
            "input": r"", # smiles
            "positive_response": "Yes, it can activate the Androgen Receptor.",
            "negative_response": "No, it cannot activate the Androgen Receptor.",
        },
    },
    "NR-AR-LBD": r"The compound {smiles} {label} activate the Androgen Receptor Ligand Binding Domain(AR-LBD).",
    "NR-AhR": r"The compound {smiles} {label} activate the aryl hydrocarbon Receptor(AhR).",
    "NR-Aromatase": r"The compound {smiles} {label} activate the aromatase.",
    "NR-ER": r"The compound {smiles} {label} activate the Estrogen Receptor(ER).",
    "NR-ER-LBD": r"The compound {smiles} {label} activate the Estrogen Receptor Ligand Binding Domain(ER-LBD).",
    "NR-PPAR-gamma": r"The compound {smiles} {label} activate the peroxisome proliferator-activated receptor gamma(PPAR-gamma).",
    "SR-ARE": r"The compound {smiles} {label} activate the Antioxidant Response Element(ARE).",
    "SR-ATAD5": r"The compound {smiles} {label} activate the ATPase family AAA domain-containing protein 5(ATAD5).",
    "SR-HSE": r"The compound {smiles} {label} activate the Heat Shock Response Element(HSE).",
    "SR-MMP": r"The compound {smiles} {label} activate the Matrix Metalloproteinase(MMP).",
    "SR-p53": r"The compound {smiles} {label} activate the tumor protein p53.",
}

QM7B = {}
QM8 = {}
QM9 = {}

DRUG1_DRUG2_Y = {
    "TWOSIDES": {
        "common": r"The concurrent use of {smiles_1} and {smiles_2} can cause many polypharmacy side-effects, including {label}.",
        "instruction_version": {
            "instruction": "Can the concurrent use of the following two drugs represented by SMILES cause side effects?",
            "input": r"smiles1: {smiles_1}, smiles2: {smiles_2},", # smiles
            "positive_response": "Yes, it can cause many side effects, including {label}.",
            "negative_response": "No, you can use it without worries.",
        },
    },
    "DrugBank": { # {1: '#Drug1 may increase the photosensitizing activities of #Drug2.', 2: '#Drug1 may increase the anticholinergic activities of #Drug2.', 3: 'The bioavailability of #Drug2 can be decreased when combined with #Drug1.', 4: 'The metabolism of #Drug2 can be increased when combined with #Drug1.', 5: '#Drug1 may decrease the vasoconstricting activities of #Drug2.', 6: '#Drug1 may increase the anticoagulant activities of #Drug2.', 7: '#Drug1 may increase the ototoxic activities of #Drug2.', 8: 'The therapeutic efficacy of #Drug2 can be increased when used in combination with #Drug1.', 9: '#Drug1 may increase the hypoglycemic activities of #Drug2.', 10: '#Drug1 may increase the antihypertensive activities of #Drug2.', 11: 'The serum concentration of the active metabolites of #Drug2 can be reduced when #Drug2 is used in combination with #Drug1 resulting in a loss in efficacy.', 12: '#Drug1 may decrease the anticoagulant activities of #Drug2.', 13: 'The absorption of #Drug2 can be decreased when combined with #Drug1.', 14: '#Drug1 may decrease the bronchodilatory activities of #Drug2.', 15: '#Drug1 may increase the cardiotoxic activities of #Drug2.', 16: '#Drug1 may increase the central nervous system depressant (CNS depressant) activities of #Drug2.', 17: '#Drug1 may decrease the neuromuscular blocking activities of #Drug2.', 18: '#Drug1 can cause an increase in the absorption of #Drug2 resulting in an increased serum concentration and potentially a worsening of adverse effects.', 19: '#Drug1 may increase the vasoconstricting activities of #Drug2.', 20: '#Drug1 may increase the QTc-prolonging activities of #Drug2.', 21: '#Drug1 may increase the neuromuscular blocking activities of #Drug2.', 22: '#Drug1 may increase the adverse neuromuscular activities of #Drug2.', 23: '#Drug1 may increase the stimulatory activities of #Drug2.', 24: '#Drug1 may increase the hypocalcemic activities of #Drug2.', 25: '#Drug1 may increase the atrioventricular blocking (AV block) activities of #Drug2.', 26: '#Drug1 may decrease the antiplatelet activities of #Drug2.', 27: '#Drug1 may increase the neuroexcitatory activities of #Drug2.', 28: '#Drug1 may increase the dermatologic adverse activities of #Drug2.', 29: '#Drug1 may decrease the diuretic activities of #Drug2.', 30: '#Drug1 may increase the orthostatic hypotensive activities of #Drug2.', 31: 'The risk or severity of hypertension can be increased when #Drug2 is combined with #Drug1.', 32: '#Drug1 may increase the sedative activities of #Drug2.', 33: 'The risk or severity of QTc prolongation can be increased when #Drug1 is combined with #Drug2.', 34: '#Drug1 may increase the immunosuppressive activities of #Drug2.', 35: '#Drug1 may increase the neurotoxic activities of #Drug2.', 36: '#Drug1 may increase the antipsychotic activities of #Drug2.', 37: '#Drug1 may decrease the antihypertensive activities of #Drug2.', 38: '#Drug1 may increase the vasodilatory activities of #Drug2.', 39: '#Drug1 may increase the constipating activities of #Drug2.', 40: '#Drug1 may increase the respiratory depressant activities of #Drug2.', 41: '#Drug1 may increase the hypotensive and central nervous system depressant (CNS depressant) activities of #Drug2.', 42: 'The risk or severity of hyperkalemia can be increased when #Drug1 is combined with #Drug2.', 43: 'The protein binding of #Drug2 can be decreased when combined with #Drug1.', 44: '#Drug1 may increase the central neurotoxic activities of #Drug2.', 45: '#Drug1 may decrease effectiveness of #Drug2 as a diagnostic agent.', 46: '#Drug1 may increase the bronchoconstrictory activities of #Drug2.', 47: 'The metabolism of #Drug2 can be decreased when combined with #Drug1.', 48: '#Drug1 may increase the myopathic rhabdomyolysis activities of #Drug2.', 49: 'The risk or severity of adverse effects can be increased when #Drug1 is combined with #Drug2.', 50: 'The risk or severity of heart failure can be increased when #Drug2 is combined with #Drug1.', 51: '#Drug1 may increase the hypercalcemic activities of #Drug2.', 52: '#Drug1 may decrease the analgesic activities of #Drug2.', 53: '#Drug1 may increase the antiplatelet activities of #Drug2.', 54: '#Drug1 may increase the bradycardic activities of #Drug2.', 55: '#Drug1 may increase the hyponatremic activities of #Drug2.', 56: 'The risk or severity of hypotension can be increased when #Drug1 is combined with #Drug2.', 57: '#Drug1 may increase the nephrotoxic activities of #Drug2.', 58: '#Drug1 may decrease the cardiotoxic activities of #Drug2.', 59: '#Drug1 may increase the ulcerogenic activities of #Drug2.', 60: '#Drug1 may increase the hypotensive activities of #Drug2.', 61: '#Drug1 may decrease the stimulatory activities of #Drug2.', 62: 'The bioavailability of #Drug2 can be increased when combined with #Drug1.', 63: '#Drug1 may increase the myelosuppressive activities of #Drug2.', 64: '#Drug1 may increase the serotonergic activities of #Drug2.', 65: '#Drug1 may increase the excretion rate of #Drug2 which could result in a lower serum level and potentially a reduction in efficacy.', 66: 'The risk or severity of bleeding can be increased when #Drug1 is combined with #Drug2.', 67: '#Drug1 can cause a decrease in the absorption of #Drug2 resulting in a reduced serum concentration and potentially a decrease in efficacy.', 68: '#Drug1 may increase the hyperkalemic activities of #Drug2.', 69: '#Drug1 may increase the analgesic activities of #Drug2.', 70: 'The therapeutic efficacy of #Drug2 can be decreased when used in combination with #Drug1.', 71: '#Drug1 may increase the hypertensive activities of #Drug2.', 72: '#Drug1 may decrease the excretion rate of #Drug2 which could result in a higher serum level.', 73: 'The serum concentration of #Drug2 can be increased when it is combined with #Drug1.', 74: '#Drug1 may increase the fluid retaining activities of #Drug2.', 75: 'The serum concentration of #Drug2 can be decreased when it is combined with #Drug1.', 76: '#Drug1 may decrease the sedative activities of #Drug2.', 77: 'The serum concentration of the active metabolites of #Drug2 can be increased when #Drug2 is used in combination with #Drug1.', 78: '#Drug1 may increase the hyperglycemic activities of #Drug2.', 79: '#Drug1 may increase the central nervous system depressant (CNS depressant) and hypertensive activities of #Drug2.', 80: '#Drug1 may increase the hepatotoxic activities of #Drug2.', 81: '#Drug1 may increase the thrombogenic activities of #Drug2.', 82: '#Drug1 may increase the arrhythmogenic activities of #Drug2.', 83: '#Drug1 may increase the hypokalemic activities of #Drug2.', 84: '#Drug1 may increase the vasopressor activities of #Drug2.', 85: '#Drug1 may increase the tachycardic activities of #Drug2.', 86: 'The risk of a hypersensitivity reaction to #Drug2 is increased when it is combined with #Drug1.'}
        "common": r"xxxxxxxxxxx", # label_list
        "instruction_version": {
            "instruction": "Is there any interaction between the following two drugs represented by SMILES?",
            "input": r"smiles1: {smiles_1}, smiles2: {smiles_2},", # smiles
            "positive_response": "",
            "negative_response": "",
        },
    },
} 

P_R = {
    "USPTO-50K": r"To synthesize {product}, we can use the reactant {reactant}.",
    "USPTO": {
        "common": r"To synthesize {product}, we can use the reactant {reactant}.",
        "instruction_version": {
            "instruction": "What can I use as reactants to synthesize the compound represented by the SMILES below?",
            "input": r"{smiles}",
            "positive_response": "Yes, it can cause many side effects, including {label}.",
            "negative_response": "No, you can use it without worries",
        },
    }
}

# R_C_P_L = {
#     "Buchwald-Hartwig": r"In Pd-catalysed Buchwald-Hartwig C-N cross coupling reactions, the yields of product {product},", # todo
#     "USPTO_Yields": r"The yields of product {product} {reactant} {catalyst} is {label}.", # todo
# }

# R_C_P = {
#     "USPTO_Catalyst": r"",
# }

# R_P = {
#     "USPTO_reaction": r"", # noticed
# }

SINGLE_REGRESSION_TASK = ["Caco2_Wang", "Lipophilicity_AstraZeneca", "Solubility_AqSolDB",
                          "HydrationFreeEnergy_FreeSolv", "PPBR_AZ", "VDss_Lombardo", "Half_Life_Obach", "Clearance_Hepatocyte_AZ",
                          # ------ ADME
                          "LD50_Zhu", "hERG_at_1uM", "hERG_at_10uM",
                          # ------ Tox
                          # ------ HTS(none)
                          # ------ QM
                          # ------ Yields
                          ]
SINGLE_CLASSIFICATION_TASK = ["PAMPA_NCATS", "HIA_Hou", "Pgp_Broccatelli", "Bioavailability_Ma", "BBB_Martins",
                              "CYP2C19_Veith", "CYP2D6_Veith", "CYP3A4_Veith", "CYP1A2_Veith", "CYP2C9_Veith",
                              "CYP2C9_Substrate_CarbonMangels", "CYP2D6_Substrate_CarbonMangels", "CYP3A4_Substrate_CarbonMangels",
                              # ------ ADME
                              "hERG", "hERG_inhib", "hERG_Karim", "AMES", "DILI", "Skin Reaction", "Carcinogens_Lagunin", "ClinTox",
                              "NR-AR", "NR-AR-LBD", "NR-AhR", "NR-Aromatase", "NR-ER", "NR-ER-LBD", "NR-PPAR-gamma",
                              "SR-ARE", "SR-ATAD5", "SR-HSE", "SR-MMP", "SR-p53",
                              # ------ Tox
                              "SARSCoV2_Vitro_Touret", "SARSCoV2_3CLPro_Diamond", "HIV",
                              "orexin1_receptor_butkiewicz", "m1_muscarinic_receptor_agonists_butkiewicz", "m1_muscarinic_receptor_antagonists_butkiewicz",
                              "potassium_ion_channel_kir2", "kcnq2_potassium_channel_butkiewicz", "cav3_t-type_calcium_channels_butkiewicz",
                              "choline_transporter_butkiewicz", "serine_threonine_kinase_33_butkiewicz", "tyrosyl-dna_phosphodiesterase_butkiewicz",
                              # ------ HTS
                              # # ------ QM(none)
                              # ------ Yields(none)
                              ]

CLASSIFICATION_LABEL = {
    "PAMPA_NCATS": {0: "low-to-moderate", 1: "high"},
    "HIA_Hou": {0: "inactive", 1: "active"},
    "Pgp_Broccatelli": {0: "inactive", 1: "active"},
    "Bioavailability_Ma": {0: "inactive", 1: "active"},
    "BBB_Martins": {0: "inactive", 1: "active"},
    "CYP2C19_Veith": {0: "inactive", 1: "active"},
    "CYP2D6_Veith": {0: "inactive", 1: "active"},
    "CYP3A4_Veith": {0: "inactive", 1: "active"},
    "CYP1A2_Veith": {0: "inactive", 1: "active"},
    "CYP2C9_Veith": {0: "inactive", 1: "active"},
    "CYP2C9_Substrate_CarbonMangels": {0: "is not", 1: "is"},
    "CYP2D6_Substrate_CarbonMangels": {0: "is not", 1: "is"},
    "CYP3A4_Substrate_CarbonMangels": {0: "is not", 1: "is"}, #todo 0/1 =>'is/is not' or 'is not/is'
    # ------ ADME
    "hERG": {0: "does not block", 1: "blocks"},
    "hERG_Karim": {0: "does not block", 1: "blocks"},
    "hERG_inhib": {0: "does not block", 1: "blocks"},
    "AMES": {0: "is not", 1: "is"},
    "DILI": {0: "cannot", 1: "can"},
    "Skin Reaction": {0: "cannot", 1: "can"},
    "Carcinogens_Lagunin": {0: "cannot", 1: "can"},
    "ClinTox": {0: "not failed", 1: "failed"},
    "NR-AR": {0: "cannot", 1: "can"},
    "NR-AR-LBD": {0: "cannot", 1: "can"},
    "NR-AhR": {0: "cannot", 1: "can"},
    "NR-Aromatase": {0: "cannot", 1: "can"},
    "NR-ER": {0: "cannot", 1: "can"},
    "NR-ER-LBD": {0: "cannot", 1: "can"},
    "NR-PPAR-gamma": {0: "cannot", 1: "can"},
    "SR-ARE": {0: "cannot", 1: "can"},
    "SR-ATAD5": {0: "cannot", 1: "can"},
    "SR-HSE": {0: "cannot", 1: "can"},
    "SR-MMP": {0: "cannot", 1: "can"},
    "SR-p53": {0: "cannot", 1: "can"},
    # ------ Tox
    "SARSCoV2_Vitro_Touret": {0: "inactive", 1: "active"},
    "SARSCoV2_3CLPro_Diamond": {0: "inactive", 1: "active"},
    "HIV": {0: "inactive", 1: "active"},
    "orexin1_receptor_butkiewicz": {0: "inactive", 1: "active"},
    "m1_muscarinic_receptor_agonists_butkiewicz": {0: "inactive", 1: "active"},
    "m1_muscarinic_receptor_antagonists_butkiewicz": {0: "inactive", 1: "active"},
    "potassium_ion_channel_kir2": {0: "inactive", 1: "active"},
    "kcnq2_potassium_channel_butkiewicz": {0: "inactive", 1: "active"},
    "cav3_t-type_calcium_channels_butkiewicz": {0: "inactive", 1: "active"},
    "choline_transporter_butkiewicz": {0: "inactive", 1: "active"},
    "serine_threonine_kinase_33_butkiewicz": {0: "inactive", 1: "active"},
    "tyrosyl-dna_phosphodiesterase_butkiewicz": {0: "inactive", 1: "active"},
    # ------ HTS
}


def drug_preprocess(drug):
    '''
    list of smiles
    '''
    return drug

def smi_preprocess(smi):
    '''
    single smiles
    '''
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return None
    
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    
    smi = Chem.MolToSmiles(mol)
    if smi is None:
        return None

    return smi

def task_hub(dataset, data, subtask=None, label_index=None):
    '''
    subtask: for multi-subtask tasks, e.g. tox21.
    label_index: for DDI task
    '''
    drug_y_task = [k for k, _ in DRUG_Y.items()]

    if (dataset in drug_y_task) and (dataset in SINGLE_REGRESSION_TASK):
        drug, y = data['Drug'], data['Y']
        prompt = DRUG_Y[dataset]

        output = get_regression_prompt(prompt, drug, y)

    elif (dataset in drug_y_task) and (dataset in SINGLE_CLASSIFICATION_TASK):
        drug, y = data['Drug'], data['Y']
        prompt = DRUG_Y[dataset]
        label_dict = CLASSIFICATION_LABEL[dataset]

        output = get_classification_prompt(prompt, drug, y, label_dict=label_dict)

    elif dataset == "herg_central":
        if subtask in HERG_CENTRAL.keys():
            drug, y = data['Drug'], data['Y']
            prompt = HERG_CENTRAL[subtask]

            if subtask == "hERG_inhib": # classification
                label_dict = CLASSIFICATION_LABEL[subtask]
                output = get_classification_prompt(prompt, drug, y, label_dict=label_dict)
            else: # regression
                output = get_regression_prompt(prompt, drug, y)
        
        else:
            print(f"Subtask \"{subtask}\" in dataset \"{dataset}\" not exist!")
            return None
        
    elif dataset == "Tox21":
        if subtask in TOX21.keys():
            drug, y = data['Drug'], data['Y']
            prompt = TOX21[subtask]
            label_dict = CLASSIFICATION_LABEL[subtask]

            output = get_classification_prompt(prompt, drug, y, label_dict=label_dict)

        else:
            print(f"Subtask \"{subtask}\" in dataset \"{dataset}\" not exist!")
            return None

    elif dataset == "DrugBank":
        drug1, drug2, y = data['Drug1'], data['Drug2'], data['Y']
        prompt_dict = label_index
        output = get_drugbank_prompt(drug1, drug2, y, prompt_dict)

    elif dataset == "TWOSIDES":
        drug1, drug2, y = data['Drug1'], data['Drug2'], data['Y']
        side_effect_dict = label_index
        prompt = DRUG1_DRUG2_Y[dataset]
        output = get_twosides_prompt(prompt, drug1, drug2, y, side_effect_dict)

    elif dataset in P_R.keys():
        product, reactant = data['input'], data['output']
        prompt = P_R[dataset]
        output = get_P_R_prompt(prompt, product, reactant)

    else:
        print(f"Dataset \"{dataset}\" not exist!")
        return None
    
    return output


def get_regression_prompt(prompt, Drug, Y):
    '''
    Regression Task
    '''
    assert len(Drug) == len(Y)

    Drug = drug_preprocess(Drug)

    output = []
    for x, y in zip(Drug, Y):
        x = smi_preprocess(x)
        if x == None:
            continue

        data_dict = {}
        data_dict["text"] = prompt.format(smiles="<<|mol0|>>", label=y)
        data_dict["entities"] = {"<<|mol0|>>": {"smiles": x}}
        output.append(json.dumps(data_dict))

    # output = [prompt.format(smiles=x, label=y) for x, y in zip(Drug, Y)]
    return output

def get_classification_prompt(prompt, Drug, Y, label_dict=None):
    '''
    Classification Task
    '''
    assert len(Drug) == len(Y)
    
    Drug = drug_preprocess(Drug)
    Y = [int(y) for y in Y]

    output = []
    for x, y in zip(Drug, Y):
        x = smi_preprocess(x)
        if x == None:
            continue
        
        data_dict = {}
        data_dict["text"] = prompt.format(smiles="<<|mol0|>>", label=label_dict[y])
        data_dict["entities"] = {"<<|mol0|>>": {"smiles": x}}
        output.append(json.dumps(data_dict))

    # output = [prompt.format(smiles=x, label=label_dict[y]) for x, y in zip(Drug, Y)]
    return output

def get_drugbank_prompt(Drug1, Drug2, Y, prompt_dict):
    assert len(Drug1) == len(Y)
    assert len(Drug2) == len(Y)

    Drug1, Drug2 = drug_preprocess(Drug1), drug_preprocess(Drug2)
    Y = [int(y) for y in Y]

    output = []
    for d1, d2, y in zip(Drug1, Drug2, Y):
        d1, d2 = smi_preprocess(d1), smi_preprocess(d2)
        if d1 == None or d2 == None:
            continue

        data_dict = {}
        prompt = prompt_dict[y]
        data_dict["text"] = prompt.replace("#Drug1", "<<|mol0|>>").replace("#Drug2", "<<|mol1|>>")
        data_dict["entities"] = {"<<|mol0|>>": {"smiles": d1}, "<<|mol1|>>": {"smiles": d2}}
        output.append(json.dumps(data_dict))
    
    return output

def get_twosides_prompt(prompt, Drug1, Drug2, Y, side_effect_dict):
    assert len(Drug1) == len(Y)
    assert len(Drug2) == len(Y)

    Drug1, Drug2 = drug_preprocess(Drug1), drug_preprocess(Drug2)
    Y = [int(y) for y in Y]

    drugpair_effect = {}
    for d1, d2, y in zip(Drug1, Drug2, Y):
        drug_pair1 = (d1, d2)
        drug_pair2 = (d2, d1)
        if drug_pair1 in drugpair_effect.keys():
            drugpair_effect[drug_pair1].append(y)
        elif drug_pair2 in drugpair_effect.keys():
            drugpair_effect[drug_pair2].append(y)
        else:
            drugpair_effect[drug_pair1] = [y]
    # print(len(drugpair_effect))
    
    output = []
    for drug_pair, y in tqdm(drugpair_effect.items()):
        d1, d2 = drug_pair[0], drug_pair[1]
        d1, d2 = smi_preprocess(d1), smi_preprocess(d2)
        if d1 == None or d2 == None:
            continue
        
        if len(y) == 1:
            side_effects = side_effect_dict[y[0]]
        else:
            side_effects = ", ".join([side_effect_dict[yy] for yy in y[:-1]]) + " and " + side_effect_dict[y[-1]]
        # side_effects = ", ".join([side_effect_dict[yy] for yy in y])

        data_dict = {}
        data_dict["text"] = prompt.format(smiles_1="<<|mol0|>>", smiles_2="<<|mol1|>>", label=side_effects)
        data_dict["entities"] = {"<<|mol0|>>": {"smiles": d1}, "<<|mol1|>>": {"smiles": d2}}
        output.append(json.dumps(data_dict))
    
    # output = [prompt.format(smiles_1=d1, smiles_2=d2, label=side_effect_dict[y]) for d1, d2, y in zip(Drug1, Drug2, Y)]
    return output

def get_P_R_prompt(prompt, P, R):
    assert len(P) == len(R)

    P, R = drug_preprocess(P), drug_preprocess(R)

    output = []
    for p, r in zip(P, R):
        p, r = smi_preprocess(p), smi_preprocess(r)
        if p == None or r == None:
            continue
        data_dict = {}
        data_dict["text"] = prompt.format(product="<<|mol0|>>", reactant=r)
        data_dict["entities"] = {"<<|mol0|>>": {"smiles": p}}

        # data_dict["text"] = prompt.format(product=p, reactant="<<|mol0|>>")
        # data_dict["entities"] = {"<<|mol0|>>": {"smiles": r}}

        output.append(json.dumps(data_dict))

    # output = [prompt.format(product=p, reactant=r) for p, r in zip(P, R)]
    return output

if __name__ == '__main__':
    print(f"Total single regression task: {len(SINGLE_REGRESSION_TASK)}")
    print(f"Total single classification task: {len(SINGLE_CLASSIFICATION_TASK)}")