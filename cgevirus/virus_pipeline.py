import os
import sys
import subprocess
import csv
import gzip

from cgevirus import kma

def virus_pipeline(args):
    print("Starting the virus pipeline...")

    # Check if output folder already exists
    output_dir = '/var/lib/cge/results/{}'.format(args.name)
    if os.path.exists(output_dir):
        sys.exit(
            f"Error: Output directory '{output_dir}' already exists. Please choose a different name or delete the existing directory.")

    if args.db_dir is None:
        if not os.path.exists('/var/lib/cge/database/cge_db'):
            sys.exit('Please install the cge_db. It should be located in /var/lib/cge/database/cge_db')
        else:
            args.db_dir = '/var/lib/cge/database/cge_db'
            print(f"Using CGE database directory: {args.db_dir}")

    if args.output is None:
        args.output = '/var/lib/cge/results/{}'.format(args.name)

    print(f"Creating output directory: {args.output}")
    os.system('mkdir -p ' + args.output)

    print(f"Running KMA for virus alignment on input: {args.input}")
    kma.KMARunner(args.input,
                  args.output + "/virus_alignment",
                  args.db_dir + '/virus_db/virus_db',
                  "-ont -ca -1t1 -mem_mode -t 8").run()

    print("Running KMA for CDD...")
    kma.KMARunner(args.input,
                  args.output + "/cdd",
                  args.db_dir + '/cdd_db/cdd_db',
                  "-ont -ca -1t1 -mem_mode -t 8").run()

    print("Running Prokka for annotation...")
    cmd = 'prokka -outdir {}/ --centre virus_alignment --kingdom Viruses --prefix prokka_results {}/virus_alignment.fsa --force'.format(args.output, args.output)
    os.system(cmd)

    print("Creating virus pipeline report...")
    report = create_virus_report(args)
    with open(args.output + "/report.txt", "w") as f:
        f.write(report)

    print("Virus pipeline completed successfully. Report generated and stored at {}.".format(args.output + "/report.txt"))
    return 'virus_pipeline'


def merge_fastq_files_unix(source_directory, output_name):
    """
    Merge all fastq.gz files in the given directory using Unix commands and save the output with the specified name in the home directory.

    Args:
    source_directory (str): Path to the directory containing fastq.gz files.
    output_name (str): Name for the output file.
    """
    # Home directory path
    home_directory = os.path.expanduser('~')

    # Output file path with the specified name
    output_file = os.path.join(home_directory, f'{output_name}.fastq.gz')

    # Creating the Unix command for concatenation
    cmd = f'cat {source_directory}/*.fastq.gz > {output_file}'

    # Executing the command
    subprocess.run(cmd, shell=True, check=True)

    print(f"All files merged into {output_file}")

def create_virus_report(args):
    cdd_results = read_tab_separated_file(args.output + "/cdd.res")

    highest_scoring_hit_details = get_highest_scoring_hit_details(args.output + "/virus_alignment.res")

    # Updating the report section
    report = "Virus Pipeline Results Report\n"
    report += "=" * 60 + "\n"

    # Top Scoring Virus Alignment Hit
    report += "Top Scoring Virus Alignment Hit:\n"
    report += "-" * 60 + "\n"
    if highest_scoring_hit_details:
        report += f"Template: {highest_scoring_hit_details['#Template']}\n"
        report += f"Identity: {highest_scoring_hit_details['Template_Identity']}, "
        report += f"Coverage: {highest_scoring_hit_details['Template_Coverage']}, "
        report += f"Depth: {highest_scoring_hit_details['Depth']}\n\n"
    else:
        report += "No virus alignment hits found.\n\n"

    # Prokka Results Section
    prokka_file = args.output + "/prokka_results.tsv"
    if os.path.exists(prokka_file):
        report += format_prokka_results(prokka_file)
    else:
        report += "Prokka results not found.\n"

    # CDD Results Section
    report += format_results_section(cdd_results, "CDD Findings")

    return report

def format_prokka_results(prokka_file):
    prokka_results = read_tab_separated_file(prokka_file)  # Removed the delimiter argument
    report = "Prokka Results:\n"
    report += "-" * 60 + "\n"
    for result in prokka_results:
        report += f"Locus Tag: {result['locus_tag']}, Type: {result['ftype']}, Length: {result['length_bp']} bp, Gene: {result.get('gene', 'N/A')}, EC Number: {result.get('EC_number', 'N/A')}, COG: {result.get('COG', 'N/A')}, Product: {result['product']}\n"
    report += "\n"
    return report

def read_tab_separated_file(file_path):
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        return list(reader)

def format_results_section(results, section_title):
    report = f"{section_title}:\n"
    report += "-" * 60 + "\n"
    for result in results:
        report += f"Template: {result['#Template']}\n"
        report += f"Identity: {result['Template_Identity'].strip()}, Coverage: {result['Template_Coverage'].strip()}, Depth: {result['Depth'].strip()}\n\n"
    return report

def get_highest_scoring_hit_template(file_path):
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        highest_scoring_hit = None
        max_score = float('-inf')  # Initialize with the smallest possible number

        for row in reader:
            try:
                score = float(row['Score'])
                if score > max_score:
                    highest_scoring_hit = row
                    max_score = score
            except ValueError:
                # Handle the case where the score is not a number
                continue

    return highest_scoring_hit['#Template'] if highest_scoring_hit else None


def get_highest_scoring_hit_details(file_path):
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        highest_scoring_hit = None
        max_score = float('-inf')

        for row in reader:
            try:
                score = float(row['Score'])
                if score > max_score:
                    highest_scoring_hit = row
                    max_score = score
            except ValueError:
                continue

    return highest_scoring_hit