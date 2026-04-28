import itertools
import subprocess
from concurrent.futures import ThreadPoolExecutor
import os
from tqdm import tqdm

# Load isolate names from file
ISOLATES = [
    i for i in open("rep_isolates.txt").read().split('\n') if i
]
# Define read counts and proportions
readcnts = [200]
proportions = [(0,1),(0.1,0.9),(0.2, 0.8),(0.3,0.7),(0.4,0.6),(0.5, 0.5)]
# Primer file
primer_bed_file = "primer_v3_400.bed"
ref = "reference.fasta"

# Generate isolate combinations
isolate_combinations = list(itertools.combinations(ISOLATES, 2))

# Function to execute the command for a specific set of parameters
def run_simulation(params):
    readcnt, isolate_combination, proportion = params
    proportion1, proportion2 = proportion
    file1_path = f"assemblies/{isolate_combination[0]}/{isolate_combination[0]}.fna"
    file2_path = f"assemblies/{isolate_combination[1]}/{isolate_combination[1]}.fna"
    output_path = (
        f"simulations/combinations/{readcnt}/"
        f"{isolate_combination[0]}_{proportion1}_{isolate_combination[1]}_{proportion2}"
    )
    reads_file_path = os.path.join(output_path, "reads_1.fastq")
    
    # Check if the output folder exists and contains the reads_1.fastq file
    if os.path.exists(reads_file_path):
        return "skipped"
    
    # Command to run the simulation
    command = [
        "bygul", "simulate-proportions",
        f"{file1_path},{file2_path}",
        primer_bed_file,
        ref,
        "--proportions", f"{proportion1},{proportion2}",
        "--readcnt", str(readcnt),
        "--outdir", output_path,
        "--simulator", "mason"
    ]
    try:
        subprocess.run(command, check=True)
        return "completed"
    except subprocess.CalledProcessError as e:
        print(f"Error in simulation for {isolate_combination}: {e}")
        return "error"

# Create a list of all parameter combinations
all_params = [
    (readcnt, isolate_combination, proportion)
    for readcnt in readcnts
    for isolate_combination in isolate_combinations
    for proportion in proportions
]

# Track progress using tqdm
progress_bar = tqdm(total=len(all_params), desc="Simulation Progress", unit="task")

# Wrapper function to update progress bar
def run_with_progress(params):
    result = run_simulation(params)
    progress_bar.update(1)  # Update progress bar
    return result

# Run simulations in parallel using ThreadPoolExecutor
max_workers = 20  # Adjust this to the number of CPU cores you want to use
with ThreadPoolExecutor(max_workers=max_workers) as executor:
    results = list(executor.map(run_with_progress, all_params))

# Close the progress bar
progress_bar.close()

# Summary
completed = results.count("completed")
skipped = results.count("skipped")
errors = results.count("error")
print(f"\nSummary:")
print(f"  Completed: {completed}")
print(f"  Skipped: {skipped}")
print(f"  Errors: {errors}")
