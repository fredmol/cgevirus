import os
import sys
import subprocess
import csv
import gzip

from cgevirus import kma


def virus_pipeline(args):
    if args.folder is not None:
        if args.name is None:
            sys.exit('Please provide a name for the merged file')
        else:
            merge_fastq_files(args.folder)
            args.input = os.path.join(os.path.expanduser('~'), args.name)
            #args.output = os.path.join(os.path.expanduser('~'), args.name)

    if args.db_dir is None:
        if not os.path.exists('/opt/cge/cge_db'):
            sys.exit('Please install the cge_db. It should be located in /opt/cge/cge_db')
        else:
            args.db_dir = '/opt/cge/cge_db'

    os.system('mkdir ' + args.output)
    kma.KMARunner(args.input,
              args.output + "/virus_alignment",
              args.db_dir + '/virus_db/virus_db',
              "-ont -ca -1t1 -mem_mode -t 8").run()

    kma.KMARunner(args.input,
              args.output + "/cdd",
              args.db_dir + '/cdd_db/cdd_db',
              "-ont -ca -1t1 -mem_mode -t 8").run()

    cmd = 'prokka -outdir {}/ --centre virus_alignment --kingdom Viruses --prefix prokka_results {}/virus_alignment.fsa --force'.format(args.output, args.output)
    os.system(cmd)

    report = create_virus_report(args)
    with open(args.output + "/virus_pipeline_report.txt", "w") as f:
        f.write(report)

    return 'virus_pipeline'

def merge_fastq_files(source_directory, output_name):
    """
    Merge all fastq.gz files in the given directory and save the output with the specified name in the home directory.

    Args:
    source_directory (str): Path to the directory containing fastq.gz files.
    output_name (str): Name for the output file.
    """
    # Home directory path
    home_directory = os.path.expanduser('~')

    # Output file path with the specified name
    output_file = os.path.join(home_directory, f'{output_name}.fastq.gz')

    # Get a list of all fastq.gz files in the source directory
    fastq_files = [f for f in os.listdir(source_directory) if f.endswith('.fastq.gz')]

    # Open the output file in write mode
    with gzip.open(output_file, 'wb') as f_out:
        # Iterate over each file and append its content to the output file
        for file in fastq_files:
            file_path = os.path.join(source_directory, file)
            with gzip.open(file_path, 'rb') as f_in:
                # Copy the content of each file to the output file
                f_out.writelines(f_in)

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