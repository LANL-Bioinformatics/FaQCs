#!/usr/bin/env python3
__version__ = "0.1.0"
__author__    = "Chienchi Lo, Bioscience Division, Los Alamos National Laboratory"
__date__      = "2019/08/09"
__license__ = "GPLv3"

import sys, os, errno, argparse, subprocess, shutil, datetime, time
import gzip
import logging
import log
import pandas as pd
import ugetk

script_path = os.path.dirname(os.path.abspath(__file__))

class SmartFormatter(argparse.HelpFormatter):
    def _split_lines(self, text, width):
	    if text.startswith('R|'):
	        return text[2:].splitlines()  
	    # this is the RawTextHelpFormatter._split_lines
	    return argparse.HelpFormatter._split_lines(self, text, width)

def setup_argparse():
    parser = argparse.ArgumentParser(prog='batch_FaQCs.py',
        description = '''Script to run FaQCs in batch mode''',
        formatter_class = SmartFormatter)

    inGrp = parser.add_argument_group('Required Input')
    inGrp.add_argument('-m', '--mappingFile',required=True, type=str,help='MAPPING file.  #SampleID Files.  PE reads separated by comma')   
    inGrp.add_argument('-d', '--dir', metavar='[PATH]',type=str,help='''Contains fastq files. To use this option.
                            To use this option,the mapping file need a 'Files' column with filenames for each sampleID.''')
        
    outGrp = parser.add_argument_group('Output')
    outGrp.add_argument('-o', '--outdir', metavar='[PATH]',type=str, default='.', help='Output directory')
    #outGrp.add_argument('--zip',action='store_true', help='Zip output files for visualization')

    optGrp = parser.add_argument_group('Parameters')
    optGrp.add_argument('-q', '--quality', metavar='<INT>',default=5, type=int, help = 'Targets # as quality level (default 5) for trimming')
    # optGrp.add_argument('--singleOrPairedQC', metavar='<STR>', default='single', type=str, help = 'Choose whether to use single or paired QC. single or paired [defaule: single]')
    optGrp.add_argument('--5end', metavar='<INT>', default=0, type=int, help = 'Cut # bp from 5 end before quality trimming/filtering [default: 0]')
    optGrp.add_argument('--3end', metavar='<INT>', default=0, type=int, help = 'Cut # bp from 3 end before quality trimming/filtering [default: 0]')        
    optGrp.add_argument('--min_L', metavar='<INT>', default=50, type=int, help = 'Trimmed read should have to be at least this minimum length [default: 50]')    
    optGrp.add_argument('--avg_q', metavar='<INT>', default=0, type=str, help='Average quality cutoff [default: 0] no filter')
    optGrp.add_argument('-n', '--numAmbiguity', metavar='<INT>',default=2, type=int, help = 'Trimmed read has greater than or equal to this number of continuous base "N" will be discarded. [default:2, "NN"]')
    optGrp.add_argument('--lc', metavar='<FLOAT>', default=0.85, type=float, help='Low complexity filter ratio, Maximum fraction of mono-/di-nucleotide sequence. [default: 0.85]')
    optGrp.add_argument('--adapter', action='store_true', help='Trim reads with illumina adapter/primers')
    optGrp.add_argument('--discard', action='store_true', help='Output discarded reads to prefix.discard.fastq.')
    optGrp.add_argument('--trimOnly', action='store_true', help='No quality report. Output trimmed reads only.')
    optGrp.add_argument('--5trim_off', action='store_true', help='Turn off trimming from 5\'end.')
    optGrp.add_argument('--qcOnly', action='store_true', help='no Filters, no Trimming, report numbers.')
    optGrp.add_argument('--kmerRarefaction', action='store_true', help='Turn on the kmer calculation. Turn on will slow down ~10 times.')

    optGrp.add_argument('-c', '--cpus', metavar='<INT>', default=8, type=int, help='Number of CPUS')

    parser.add_argument('--version', action='version', version='%(prog)s v{version}'.format(version=__version__))
    parser.add_argument('--quiet', action='store_true', help='Keep messages in terminal minimal')
    parser.add_argument('--verbose', action='store_true', help='Show more infomration in log')
    parser.add_argument('--force', action='store_true', help='overwrite existing QC result')
    parser.add_argument('--uge', action='store_true', help='use uge for parallel QC run')
    return parser.parse_args()

def mkdir_p(directory_name):
    try:
        os.makedirs(directory_name)
    except OSError as exc: 
        if exc.errno == errno.EEXIST and os.path.isdir(directory_name):
            pass

def dependency_check(cmd):
    path = shutil.which(cmd)
    if not path:
	    sys.exit("[ERROR] Executable %s not found." % cmd)
    return path 
	
def tool_version_check(cmd):
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    outs, errs = proc.communicate()
    return outs.decode().rstrip()

def get_runtime(start):
    end = time.time()
    hours, rem = divmod(float(end-start), 3600)
    minutes, seconds = divmod(rem, 60)
    runtime ="Running time: {:0>2}:{:0>2}:{:0>2}".format(int(hours),int(minutes), int(seconds))
    return runtime

def process_cmd(cmd, msg=''):
    if msg:
        logger.info('FaQCs on sample ' + msg)

    logger.debug("FaQCs CMD: %s" %(cmd) )
    start=time.time()

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    outs, errs = proc.communicate()
    if proc.returncode != 0: 
        logger.error("Failed %d %s %s" % (proc.returncode, outs, errs))

    if msg:
        logger.info("FaQCs %s" % get_runtime(start))  
    
    return 0
            
def print_parameters(argvs,FaQCs_path):
    logger.info("Batch FaQCs v%s Start", __version__)
    arguments_info="""Arguments:
        Mapping File    :  %s
        Input Path      :  %s
        Output Path     :  %s    
        Quality Level   :  %d
        5'end trim      :  %d bp
        3'end trim      :  %d bp
        Minimum Len     :  %d
        Mean quality    :  %d
        Num Ambiguity   :  %d
        Low Complexity  :  %.2f
        Filter Adapter  :  %s
        Output discard  :  %s
        Trim only       :  %s
        5'Trim off      :  %s
        QC only         :  %s
        Kmer Rarefaction:  %s
        CPU Number      :  %d
        FaQCs Path      :  %s
        FaQCs Version   :  %s
    """ % (argvs.mappingFile, 
           argvs.dir,
           os.path.abspath(argvs.outdir), 
           argvs.quality,
           getattr(argvs, '5end'),
           getattr(argvs, '3end'),
           argvs.min_L,
           argvs.avg_q,
           argvs.numAmbiguity,
           argvs.lc,
           argvs.adapter,
           argvs.discard,
           argvs.trimOnly,
           getattr(argvs, '5trim_off'),
           argvs.qcOnly,
           argvs.kmerRarefaction,
           argvs.cpus, 
           dependency_check("FaQCs"),
           tool_version_check([FaQCs_path , '--version']))

    logger.info(arguments_info)

def get_FaQCs_cmd(path, argvs):
    FaQCs_cmd = [path, '--ascii 33','-q', str(argvs.quality), '--min_L', str(argvs.min_L), '-n',str(argvs.numAmbiguity), '--lc',str(argvs.lc), '-t',str(argvs.cpus)]
    if getattr(argvs, '5end') > 0 :
        FaQCs_cmd.extend(['--5end' , str(getattr(argvs, '5end'))])
    if getattr(argvs, '3end') > 0 :
        FaQCs_cmd.extend(['--3end' , str(getattr(argvs, '3end'))])
    if argvs.avg_q > 0 :
        FaQCs_cmd.extend(['--avg_q' , str(argvs.avg_q)])
    if argvs.adapter:
        FaQCs_cmd.append('--adapter')
    if argvs.discard:
        FaQCs_cmd.append('--discard')
    if argvs.trimOnly:
        FaQCs_cmd.append('--trim_only')
    if getattr(argvs, '5trim_off'):
        FaQCs_cmd.append('--5trim_off')
    if argvs.qcOnly:
        FaQCs_cmd.append('--qc_only')
    if argvs.kmerRarefaction:
        FaQCs_cmd.append('--kmer_rarefaction')
    
    return FaQCs_cmd


def get_open_function(filename):
    """
    Returns either open or gzip.open
    """
    if filename.lower().endswith('.gz'):
        return gzip.open
    else:  # plain text
        return open
    
    
def get_sequence_file_type(filename):
        """
        Determines whether a file is FASTA or FASTQ.
        """
        open_func = get_open_function(filename)
        with open_func(filename, 'rt') as seq_file:
            try:
                first_char = seq_file.read(1)
            except UnicodeDecodeError:
                first_char = ''
    
        if first_char == '>':
            return 'FASTA'
        elif first_char == '@':
            return 'FASTQ'
    
def check_file_input(dir, filename):
    seqformat = "Unknown"
    path_to_file = os.path.join(dir,filename)
    if not os.path.isfile(path_to_file):
        logger.error('could not find ' + filename)
        
    fileformat = get_sequence_file_type(path_to_file)
    if fileformat == "FASTA" : 
        logger.warning(filename + ' is in FASTA format')
        seqformat ="FASTA"
    elif fileformat == "FASTQ" : 
        logger.debug(filename + ' is in FASTQ format')
        seqformat = "FASTQ"
    else:
        logger.error("File %s is neither FASTA or FASTQ" % filename)

    return seqformat


if __name__ == '__main__':
    argvs    = setup_argparse()
    begin_t  = time.time()

    FaQCs_Path = dependency_check("FaQCs")

    mkdir_p(argvs.outdir)
    abs_output = os.path.abspath(argvs.outdir)
    abs_input = os.path.abspath(argvs.dir)
    logger = logging.getLogger(__name__)
    log.log_init(logger,
                 os.path.join(abs_output,os.path.splitext(os.path.basename(__file__))[0]+'.log'),
                 verbose=argvs.verbose,
                 quiet=argvs.quiet)

    # arguments info
    print_parameters(argvs,FaQCs_Path)
    
    # Convert excel format to tsv 
    mappingFile = os.path.abspath(argvs.mappingFile)
    mappingFile_cp = os.path.join(abs_output,'sample_metadata.txt')
    if mappingFile.lower().endswith(('.xlsx','.xls')):
        xlsx2csv_path = dependency_check("xlsx2csv")
        process_cmd([xlsx2csv_path, '-d ','tab', mappingFile, mappingFile_cp]) 
    else:
        shutil.copy(mappingFile, mappingFile_cp) 

    os.chdir(abs_output)
    # Parse Mapping file
    df =  pd.read_csv(mappingFile_cp,sep="\t")
    if 'Files' not in df.columns:
        sys.exit("[ERROR] 'Files' column not found in meta data mapping file.")

    read1List=[]
    read2List=[]
    skip=[]
    jobids = []
    for index, row in df.iterrows():
        cmd = get_FaQCs_cmd(FaQCs_Path,argvs)
        sample_out_dir = row['#SampleID']
        out_read1 = os.path.join(abs_output,sample_out_dir,'QC.1.trimmed.fastq')
        out_read2 = os.path.join(abs_output,sample_out_dir,'QC.2.trimmed.fastq')
        out_single = os.path.join(abs_output,sample_out_dir,'QC.unpaired.trimmed.fastq')
        out_pdf = os.path.join(abs_output,sample_out_dir,'QC_qc_report.pdf')
        jid = None
        if argvs.uge:
            jobname= "FaQCs"+ row['#SampleID']
            mem='10G'
            cpu = argvs.cpus
            qlog = os.path.join(abs_output,'uge.log')
            email = None
        if (',' in row['Files']):
            f_fq,r_fq = row['Files'].split(',')
            #type='pe'
            if f_fq == r_fq:
                logger.error('first and second read pair files cannot be the same file')

            if check_file_input(abs_input, f_fq) == "FASTQ"  and check_file_input(abs_input, r_fq) == "FASTQ":
                mkdir_p(sample_out_dir)
                cmd.extend(['-d', "'{0}'".format(sample_out_dir), '-1',"'{0}'".format(os.path.join(abs_input,f_fq)),'-2',"'{0}'".format(os.path.join(abs_input,r_fq))])
                if not os.path.isfile(out_read1) or not os.path.isfile(out_read2) or not os.path.isfile(out_pdf):
                    jid = ugetk.qsub(cmd,jobname,mem,cpu,qlog,email) if argvs.uge else process_cmd(cmd,sample_out_dir)
                    #logger.info("%s Not qc" % row['#SampleID'])
                else:
                    if argvs.force :
                        jid = ugetk.qsub(cmd,jobname,mem,cpu,qlog,email) if argvs.uge else process_cmd(cmd,sample_out_dir)
                    else:
                        logger.info("Result exist. Skip FaQCs for " + sample_out_dir )

                read1List.append(out_read1)
                read2List.append(out_read2)
                skip.append(False)
                if jid is not None:
                    jobids.append(jid)
            else: # Not FASTQ pe.  Please check
                skip.append(True)
                read1List.append(os.path.join(abs_input, f_fq))
                read2List.append(os.path.join(abs_input, r_fq))
                logger.info("Not paired-end reads in FASTQ format. Skip FaQCs for " + sample_out_dir )
        else:
            f_fq = row['Files']
            read2List.append("")
            #type='se'
            #print(mem)
            se_format = check_file_input(abs_input, f_fq)
            if se_format == "FASTQ":
                mkdir_p(sample_out_dir)
                cmd.extend(['-d',"'{0}'".format(sample_out_dir),'-u',"'{0}'".format(os.path.join(abs_input,f_fq))])
                if not os.path.isfile(out_single) or not os.path.isfile(out_pdf): 
                    jid = ugetk.qsub(cmd,jobname,mem,cpu,qlog,email) if argvs.uge else process_cmd(cmd,sample_out_dir)
                    #logger.info("%s Not qc" % row['#SampleID'])
                else: 
                    if argvs.force:
                        jid = ugetk.qsub(cmd,jobname,mem,cpu,qlog,email) if argvs.uge else process_cmd(cmd,sample_out_dir)
                    else:
                        logger.info("Result exist. Skip FaQCs for " + sample_out_dir )
                read1List.append(out_single)
                skip.append(False)
                if jid is not None:
                    jobids.append(jid);
            elif se_format == "FASTA":
                skip.append(True)
                read1List.append(os.path.join(abs_input, f_fq))
                logger.info("Skip FaQCs for " + sample_out_dir )
            else:
                skip.append(True)
                read1List.append(None)
                logger.info("Skip FaQCs for " + sample_out_dir )


    logger.info("Total Num of Sample: %d" % len(df))

    if argvs.uge:
       logger.info("Wait for submitted %d jobs ..." % len(jobids))
       ugetk.qwait_all(jobids)
    df['Skip_FaQCs'] = skip
    df['read1'] = read1List
    df['read2'] = read2List
    #df.to_excel("out.xlsx")
    out_for_newfof = os.path.join(abs_output,"snippy_input.tab")
    with open (out_for_newfof, "w") as f:
        for index, row in df.iterrows():
            sample_out_dir = row['#SampleID']
            out_read1 = os.path.join(abs_output,sample_out_dir,'QC.1.trimmed.fastq')
            out_read2 = os.path.join(abs_output,sample_out_dir,'QC.2.trimmed.fastq')
            out_single = os.path.join(abs_output,sample_out_dir,'QC.unpaired.trimmed.fastq')
            out_pdf = os.path.join(abs_output,sample_out_dir,'QC_qc_report.pdf')
            #if os.path.isfile(out_read1) and os.path.isfile(out_read2) and os.path.isfile(out_pdf):
            if os.path.isfile(out_read1) and os.path.isfile(out_read2):
                w_string = "\t".join([row['#SampleID'],out_read1,out_read2])
            #elif os.path.isfile(out_read1) and os.path.isfile(out_pdf):
            elif os.path.isfile(out_single):
                w_string = "\t".join([row['#SampleID'],out_single])
            elif row['Skip_FaQCs'] is True and row['read1'] is not None:
                w_string = "\t".join([row['#SampleID'],row['read1']])
            else:
                w_string = None
                logger.warning("FaQCs may fail on %s" % row['#SampleID'])

            if w_string is not None:
                f.write(w_string + "\n")
        f.close()
        
    #df.to_csv(out_for_newfof,sep="\t",index=None, header=False,columns= ['#SampleID', 'read1','read2'])
    
    logger.info("Num of Sample by FaQCs: %d" % df.Skip_FaQCs.value_counts().loc[False])
    totol_runtime = "Total %s" % get_runtime(begin_t)
    logger.info(totol_runtime)
    if argvs.quiet:
        print(totol_runtime)
    
    
    
