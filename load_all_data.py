import argparse
import os

from load_data import (ALL_TASK, GENERATE_RETROSYN_TASK, MULTI_DDI_TASK,
                       SINGLE_TASK, get_outputs_of_dataset)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=str, default="prompt_data")
    parser.add_argument('--task-class', choices=['single', 'multi', 'generate', 'all'], default='single')
    args = parser.parse_args()

    all_task, single_task, multi_task, generate_task = ALL_TASK, SINGLE_TASK, MULTI_DDI_TASK, GENERATE_RETROSYN_TASK
    task_dict = {
        "all": all_task, "single": single_task, "multi": multi_task, "generate": generate_task,
    }

    outputs = []
    for task in task_dict[args.task_class]:
        print(f"Preprocessing dataset {task}...")
        output = get_outputs_of_dataset(task, "random")
        outputs.extend(output)
    
    # output_dir = os.path.join(args.output_dir, args.task_class)
    # if not os.path.exists(output_dir):
    #     os.makedirs(output_dir)

    with open(os.path.join(args.output_dir, f"{args.task_class}_prompts_v2.jsonl"), "w", encoding="utf8") as fw:
        for output in outputs:
            fw.write(output + "\n")

    print(f"Total outputs: {len(outputs)}.")
