import argparse
import os

from load_data import ALL_TASK, get_outputs_of_dataset

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=str, default="prompt_data")
    args = parser.parse_args()

    all_task = ALL_TASK
    outputs = []
    for task in all_task:
        print(f"Preprocessing dataset {task}...")
        output = get_outputs_of_dataset(task, "random")
        outputs.extend(output)
    
    with open(os.path.join(args.output_dir, "all_prompts_v1.txt"), "w", encoding="utf8") as fw:
        for output in outputs:
            print(output, file=fw)


