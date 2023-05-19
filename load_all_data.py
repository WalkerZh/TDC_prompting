import argparse
import os

from load_data import (ALL_TASK, CLASSIFICATION_ADME_TOX_TASK,
                       CLASSIFICATION_HTS_TASK, CLASSIFICATION_TASK,
                       GENERATE_RETROSYN_TASK, MOLECULENET_TASK,
                       MULTI_DDI_TASK, REGRESSION_TASK, SINGLE_TASK,
                       get_outputs_of_dataset)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=str, default="prompt_data")
    parser.add_argument('--task-class', choices=['single', 'multi', 'generate', 'all', 
                                                 'classification', 'regression', 'classification_ADME_TOX', 'classification_HTS'], default='single')
    args = parser.parse_args()

    all_task, single_task, multi_task, generate_task = ALL_TASK, SINGLE_TASK, MULTI_DDI_TASK, GENERATE_RETROSYN_TASK
    molnet_task, classification_task, regression_task = MOLECULENET_TASK, CLASSIFICATION_TASK, REGRESSION_TASK
    classification_adme_tox_task, classification_hts_task = CLASSIFICATION_ADME_TOX_TASK, CLASSIFICATION_HTS_TASK
    task_dict = {
        "all": all_task, "single": single_task, "multi": multi_task, "generate": generate_task,
        "classification": classification_task, "regression": regression_task,
        "classification_ADME_TOX": classification_adme_tox_task, "classification_HTS": classification_hts_task,
    }

    outputs = []
    for task in task_dict[args.task_class]:
        if task == "herg_central":
            print(f"Dataset {task} passed.")
            continue
        if task in molnet_task:
            print(f"Dataset {task} ignored because it's in MoleculeNet.")
            continue
        print(f"Preprocessing dataset {task}...")
        if args.task_class in ['classification', 'classification_ADME_TOX', 'classification_HTS']:
            output = get_outputs_of_dataset(task, "random", reg=0)
        else:
            output = get_outputs_of_dataset(task, "random")
        print(f"Total examples of dataset {task}: {len(output)}")
        outputs.extend(output)
    
    # output_dir = os.path.join(args.output_dir, args.task_class)
    # if not os.path.exists(output_dir):
    #     os.makedirs(output_dir)

    with open(os.path.join(args.output_dir, f"{args.task_class}_prompts_v3.jsonl"), "w", encoding="utf8") as fw:
        for output in outputs:
            fw.write(output + "\n")

    print(f"Total outputs: {len(outputs)}.")
