import os, sys
import glob
from csv import DictReader
from pbtranscript.Utils import real_upath

def generate_batch_cmds_for_polishing(chunk_prefix, nfl_filename, subread_xml, cpus, cmd_filename, walltime, queue):

    subread_xml = real_upath(subread_xml)
    nfl_filename = real_upath(nfl_filename)

    fastas = glob.glob(chunk_prefix + '*.consensus.fasta')
    # verify that the pickles exists as well
    for fasta in fastas:
        pickle = fasta[:-len('.consensus.fasta')] + '.pickle'
        print "looking for", pickle
        assert os.path.exists(pickle)
        dirname = fasta[:-len('.consensus.fasta')]
        if os.path.exists(dirname):
            print >> sys.stderr, "Directory {0} already exist! Abort!".format(dirname)
            sys.exit(-1)

    cmd_f = open(cmd_filename, 'w')

    for fasta in fastas:
        pickle = fasta[:-len('.consensus.fasta')] + '.pickle'
        dirname = fasta[:-len('.consensus.fasta')]
        full_fasta = real_upath(fasta)
        full_pickle = real_upath(pickle)
        os.makedirs(dirname)
        os.chdir(dirname)
        os.symlink(full_fasta, os.path.basename(full_fasta))
        os.makedirs('output')
        os.chdir('output')
        os.symlink(full_pickle, 'final.pickle')
        os.symlink(full_fasta, 'final.consensus.fasta')
        os.chdir('../../')
        f = open(os.path.join(dirname, dirname+'.sh'), 'w')
        f.write("#!/bin/bash\n")
        f.write("source /projects/banchereau-lab/ISO-seq/annotation_processing/pitchfork_ToFU2_dev/setup-env.sh\n")
        f.write("module load gcc/4.9.2\n")
        f.write("module load graphviz\n")
        f.write("PATH=$PATH:/projects/banchereau-lab/ISO-seq/annotation_processing/cDNA_Cupcake/sequence\n")
        f.write("cd $PBS_O_WORKDIR\n")
        f.write("run_IcePartial2.py all {nfl} {p}.consensus.fasta {p}.nfl.pickle "\
                "--root_dir {d} --aligner_choice=daligner --cpus={c}\n".format(\
                p=dirname, nfl=nfl_filename, d=real_upath(dirname), c=cpus))
        f.write("run_IceArrow2.py all {d} --subread_xml {s} --blasr_nproc {c} --arrow_nproc {c} --hq_min_full_length_reads=1\n".format(\
                d=real_upath(dirname), s=subread_xml, c=cpus))
        f.close()
        cmd_f.write("qsub -q {q} -l walltime={w} -l nodes=1:ppn={c} {sh}\n".format(sh=real_upath(f.name), c=cpus, w=walltime, q=queue))





if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser("Generate batch commands for running IcePartial2->IceArrow2 for each collected chunk")
    parser.add_argument("chunk_prefix", help="Collected chunk prefix (ex: collected_final.chunk)")
    parser.add_argument("nfl_filename", help="nFL read filename (ex: isoseq_nfl.fasta)")
    parser.add_argument("subread_xml", help="Subread XML file")
    parser.add_argument("--cpus", default=12, type=int, help="Number of CPUs (default: 12)")
    parser.add_argument("--cmd_filename", default='cmds', help="Output command filename (default: cmds)")
    parser.add_argument("--walltime", default='24:00:00', help="walltime for queue (default: 24:00:00)")
    parser.add_argument("--queue", default='batch2', help="queue (default: batch2)")
    args = parser.parse_args()

    generate_batch_cmds_for_polishing(args.chunk_prefix, args.nfl_filename, args.subread_xml,
                                      args.cpus, args.cmd_filename, args.walltime, args.queue)
