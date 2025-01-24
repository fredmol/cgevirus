import os
import sys
import subprocess
import csv
import gzip
import shutil

from cgevirus import kma


# pdf report
from datetime import datetime
from pathlib import Path
from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML, CSS
import base64
from cgevirus.version import __version__
import numpy as np
from io import BytesIO

import matplotlib
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['figure.figsize'] = [8, 4]
import matplotlib.pyplot as plt

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

    # Copy the input FASTQ file to the output directory
    print(f"Copying input FASTQ file to the output directory: {args.input} -> {args.output}")
    shutil.copy(args.input, args.output)

    print(f"Running KMA for virus alignment on input: {args.input}")
    kma.KMARunner(args.input,
                  args.output + "/virus_alignment",
                  args.db_dir + '/virus_db/virus_db',
                  "-ont -ca -1t1 -mem_mode -t 8 -ef").run()

    print("Running Prokka for annotation...")
    cmd = 'prokka -outdir {}/ --centre virus_alignment --kingdom Viruses --prefix prokka_results {}/virus_alignment.fsa --force'.format(args.output, args.output)
    os.system(cmd)

    print("Creating virus pipeline report...")
    report = create_virus_report(args)
    with open(args.output + "/report.txt", "w") as f:
        f.write(report)
        
    ################# PDF REPORT #################
    
    # Add these lines for PDF report generation
    print("Creating PDF report...")
    highest_scoring_hit_details = get_highest_scoring_hit_details(args.output + "/virus_alignment.res")
    
    # Read Prokka results
    prokka_results = []
    prokka_file = os.path.join(args.output, "prokka_results.tsv")
    if os.path.exists(prokka_file):
        prokka_results = read_tab_separated_file(prokka_file)

    pdf_path = create_pdf_report(args, highest_scoring_hit_details, prokka_results)
    if pdf_path:
        print(f"PDF report generated and stored at: {pdf_path}")

    ################# PDF REPORT DONE #################

    print("Virus pipeline completed successfully. Reports generated and stored at {}.".format(args.output))
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

    highest_scoring_hit_details = get_highest_scoring_hit_details(args.output + "/virus_alignment.res")

    # Updating the report section
    report = "Virus Analysis report: {}\n".format(args.name)
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

def get_virus_hits(file_path, max_hits=5):
    """Get top virus hits, sorted by score"""
    try:
        with open(file_path, 'r') as file:
            reader = csv.DictReader(file, delimiter='\t')
            hits = []
            for row in reader:
                try:
                    hits.append({
                        '#Template': row['#Template'],
                        'Score': float(row['Score']),
                        'Template_Identity': row['Template_Identity'],
                        'Template_Coverage': row['Template_Coverage'],
                        'Depth': row['Depth']
                    })
                except (ValueError, KeyError):
                    continue
            
            # Sort by score and return top hits
            hits.sort(key=lambda x: x['Score'], reverse=True)
            return hits[:max_hits]
    except Exception:
        return []


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
    
    
    
########################## PDF REPORT #########################

def get_file_stats(input_file, output_dir):
    """Get file stats using quick commands and mapstat info."""
    file_name = 'Unknown'
    file_size = 'N/A'
    fragment_count = None
    kma_version = None
    
    # Get file_name and file_size
    if os.path.exists(input_file):
        file_name = os.path.basename(input_file)
        try:
            file_size_cmd = f"ls -lh {input_file} | cut -d' ' -f5"
            file_size = subprocess.check_output(file_size_cmd, shell=True).decode().strip()
            if not file_size:
                file_size = 'N/A'
        except:
            file_size = 'N/A'
    
    # Get fragment_count and kma_version from mapstat
    mapstat_file = os.path.join(output_dir, "virus_alignment.mapstat")
    if os.path.exists(mapstat_file):
        try:
            with open(mapstat_file, 'r') as f:
                for line in f:
                    if line.startswith("## fragmentCount"):
                        parts = line.strip().split('\t')
                        if len(parts) > 1:
                            try:
                                fragment_count = int(parts[1])
                            except ValueError:
                                fragment_count = None
                    elif line.startswith("## version"):
                        parts = line.strip().split('\t')
                        if len(parts) > 1:
                            kma_version = parts[1]
        except:
            pass

    return {
        'file_name': file_name if file_name else 'Unknown',
        'file_size': file_size if file_size else 'N/A',
        'read_count': f"{fragment_count:,}" if fragment_count is not None else 'N/A',
        'kma_version': kma_version if kma_version else 'N/A'
    }

def create_pdf_report(args, highest_scoring_hit_details, prokka_results):
    """Generate a PDF report for the virus analysis"""
    try:
        # Setup paths
        package_dir = Path(__file__).parent
        assets_dir = package_dir / 'assets'
        logo_path = assets_dir / 'dtu_logo.png'

        # Prepare template data
        virus_hits = get_virus_hits(args.output + "/virus_alignment.res")
        template_data = {
            'name': getattr(args, 'name', 'Unknown'),
            'date': datetime.now().strftime("%B %d, %Y"),
            'version': __version__,
            'logo_data_url': '',
            'file_stats': get_file_stats(args.input, args.output),
            'virus_name': 'No virus identified',
            'virus_details': {
                'name': 'No virus identified',
                'identity': 'N/A',
                'coverage': 'N/A',
                'depth': 'N/A'
            },
            'identity_value': 'N/A',
            'coverage_value': 'N/A',
            'depth_value': 'N/A',
            'gene_count': len(prokka_results) if prokka_results else 0,
            'gc_content': calculate_gc_content(args.output + "/virus_alignment.fsa"),
            'prokka_results': [],
            'additional_virus_hits': [],
            'has_multiple_hits': False
        }

        if virus_hits:
            primary_hit = virus_hits[0]
            template_data['virus_name'] = primary_hit['#Template']
            template_data['identity_value'] = f"{float(primary_hit['Template_Identity']):.2f}"
            template_data['coverage_value'] = f"{float(primary_hit['Template_Coverage']):.2f}"
            template_data['depth_value'] = f"{float(primary_hit['Depth']):.2f}"
            
            if len(virus_hits) > 1:
                template_data['has_multiple_hits'] = True
                template_data['additional_virus_hits'] = [
                    {
                        'name': hit['#Template'],
                        'identity': f"{float(hit['Template_Identity']):.2f}",
                        'coverage': f"{float(hit['Template_Coverage']):.2f}",
                        'depth': f"{float(hit['Depth']):.2f}"
                    }
                    for hit in virus_hits[1:]
                ]

        
        # Load and encode logo
        if logo_path.exists():
            with open(logo_path, 'rb') as f:
                logo_data = f.read()
                template_data['logo_data_url'] = f'data:image/png;base64,{base64.b64encode(logo_data).decode("utf-8")}'

        # Format virus details
        if highest_scoring_hit_details:
            try:
                template_data['virus_details'] = {
                    'name': highest_scoring_hit_details['#Template'],
                    'identity': f"{float(highest_scoring_hit_details['Template_Identity']):.2f}",
                    'coverage': f"{float(highest_scoring_hit_details['Template_Coverage']):.2f}",
                    'depth': f"{float(highest_scoring_hit_details['Depth']):.2f}"
                }
                template_data['virus_name'] = template_data['virus_details']['name']
                template_data['identity_value'] = f"{float(highest_scoring_hit_details['Template_Identity']):.2f}"
                template_data['coverage_value'] = f"{float(highest_scoring_hit_details['Template_Coverage']):.2f}"
                template_data['depth_value'] = f"{float(highest_scoring_hit_details['Depth']):.2f}"
            except (KeyError, ValueError):
                pass

        # Format Prokka results
        if prokka_results:
            template_data['prokka_results'] = prokka_results

        # Setup Jinja2 environment and render
        env = Environment(
            loader=FileSystemLoader(package_dir / 'templates'),
            autoescape=True
        )
        template = env.get_template('report.html')
        html_content = template.render(**template_data)

        # Create PDF
        pdf_path = Path(args.output) / f"{args.name}_report.pdf"
        HTML(string=html_content).write_pdf(
            pdf_path,
            stylesheets=[CSS(assets_dir / 'style.css')]
        )

        return pdf_path

    except Exception as e:
        print(f"Error generating PDF report: {str(e)}")
        return None
        
        
def calculate_gc_content(fasta_file):
    """Calculate GC content from FASTA file"""
    try:
        if not os.path.exists(fasta_file):
            return 'N/A'
            
        gc_count = 0
        total_count = 0
        
        with open(fasta_file, 'r') as f:
            for line in f:
                if not line.startswith('>'):
                    sequence = line.strip().upper()
                    gc_count += sequence.count('G') + sequence.count('C')
                    total_count += len(sequence)
        
        if total_count > 0:
            gc_content = (gc_count / total_count) * 100
            return f"{gc_content:.1f}"
        return 'N/A'
    except:
        return 'N/A'
        