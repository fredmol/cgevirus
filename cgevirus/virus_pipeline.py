import os
import sys
import subprocess
import csv

from cgevirus import kma


def virus_pipeline(args):
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
    highest_scoring_hit = get_highest_scoring_hit_template(args.output + "/virus_alignment.res")

    kma.KMARunner(args.input,
              args.output + "/cdd",
              args.db_dir + '/cdd_db/cdd_db',
              "-ont -ca -1t1 -mem_mode -t 8").run()

    cmd = 'prokka -outdir {}/ --centre virus_alignment --kingdom Viruses --prefix prokka_results {}/virus_alignment.fsa --force'.format(args.output, args.output)
    os.system(cmd)

    create_virus_report(args, highest_scoring_hit)

    return 'virus_pipeline'

def create_virus_report(args, highest_scoring_hit):
    cdd_results = read_tab_separated_file(args.output + "/cdd.res")

    highest_scoring_hit = get_highest_scoring_hit_template(args.output + "/virus_alignment.res")

    report = "Virus Pipeline Results Report\n"
    report += "=" * 60 + "\n"

    # Top Scoring Virus Alignment Hit
    report += "Top Scoring Virus Alignment Hit:\n"
    report += "-" * 60 + "\n"
    report += f"Template: {highest_scoring_hit}\n\n"


    # Prokka Results Section
    prokka_file = os.path.join(args.output, "prokka_results", "prokka_results.tsv")
    if os.path.exists(prokka_file):
        report += format_prokka_results(prokka_file)
    else:
        report += "Prokka results not found.\n"

    # CDD Results Section
    report += format_results_section(cdd_results, "CDD Findings")

    return report

def format_prokka_results(prokka_file):
    prokka_results = read_tab_separated_file(prokka_file, delimiter='\t')
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