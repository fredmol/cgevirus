import os
import sys
import subprocess

from cgevirus import kma


def virus_pipeline(args):
    os.system('mkdir ' + args.output)
    # Check if kma is installed
    kma.KMARunner(args.input,
              args.output + "/virus_alignment",
              args.db_dir + '/virus_db/virus_db',
              "-ont -ca -1t1 -mem_mode").run()

    kma.KMARunner(args.input,
              args.output + "/cdd",
              args.db_dir + '/cdd_db/cdd_db',
              "-ont -ca -1t1 -mem_mode").run()

    cmd = 'prokka -outdir {}/ --centre virus_alignment --kingdom Viruses --prefix prokka_results {}/virus_alignment.fsa --force'.format(args.output, args.output)
    os.system(cmd)


    #pathogenicity?

    #Compile results TBD

    #Same as for cgeisolate, compile all results into a file.
        #Should be the CDD results, the prokka results, and the kma results (virus with top alignment score).


    return 'virus_pipeline'