=pod

=head1 NAME

OxfordPlots

=head1 SYNOPSIS

init_pipeline.pl OxfordPlots --pipeline_url $HIVE_URL --ref_genome homo_sapiens [options]
 
runWorker.pl -url $HIVE_URL
 
beekeeper.pl -url $HIVE_URL -loop
 
=head1 DESCRIPTION

This pipeline runs LastZ between one reference genome and all other ones. It has been used to align
the human genome to a set of other primates and to align the pig sex chromosomes to other mammalian
sex chromosomes

It assumes the following organisation of files:
 - Everything happens in a working directory call work_dir (the current directory by default
 - The input FASTA files are in a sub-directory called "genomes" (this can be configured)
 - Inside the "genomes" sub-directory, there is one directory per species. The name of that
        directory is the binomial name of the species in lowercase and with underscores instead of
        whitespaces. In addition, there is one text file per species with chromosome names and
        lengths. Only include the chromosomes you want to display. The file is named like the
        directory but ends in ".txt".
 - The reference genome (ref_genome) must be spelled like the directory name
 - All the alignments are stored in a sub-directory called "alignments" (this can be configured)

The pipeline has the following structure:

=head2 other_genomes

This is the first analysis. It reads the content of work_dir/genomes, skips the reference genome and
sends 1 job per directory (i.e. non-reference genome) to the "other_genome_files" analysis.

Runs locally.

=head2 other_genome_files

For each non-reference genome, lists the FASTA files in that directory. Send one job per FASTA file
to the "ref_genome_files" analysis.

Runs locally.

=head2 ref_genome_files

For the reference genome, lists the FASTA files in that directory. Send one job per FASTA file to
the "lastz" analysis. Because its input comes from the "other_genome_files" analysis, this creates
one "lastz" job per pair (non-reference vs reference) of chromosomes.

Runs locally.

=head2 lastz

Runs lastz with the given options (lastz_options). Stores the axt and rdotplot files in the
work_dir/alignments/non-ref-genome directory

Assumes lastz is in the path.

=head2 axt_chain

Runs axtChain with the given options (axt_chain_options) on the lastz alignments.

Assumes axtChain is in the path.

=head2 chain_to_rdotplot

Reads the chain files and produces rdotplot files, i.e. files that can be read in R for displaying
the chains on an Oxford plot for instance.

Runs locally.

=head2 make_plot

Build a dotplot (Oxfrod plot) for each pair of species. This runs an Rscript (provided with this
module) that reads the rdotplot files from the previous step and combines them in the same plot. If
you are aligning one single chromosome, this will use just one rdotplot file.

=head1 OPTIONS

=head2 pipeline_url (def: sqlite:///oxford_plots)

Defines where the eHive database lives. This can be an SQLite DB (sqlite:///db_name) or a MySQL
database (sql://user:pass@host:port/db_name)

=head2 ref_genome (mandatory)

Specifies which genome is the reference one. In other words, the pipeline will align the sequences
from all the other genomes to this genome.

=head2 work_dir (def: $CWD)

Specifies where the genomes and alignments directories are.

=head2 input_dir (def: $work_dir/genomes)

Overwrites the location of the input directory. Provide the full path.

This directory contains the fasta files organised by species. It also contains the text files with
the chromosome lengths for each species.

=head2 output_dir (def: $work_dir/alignments)

Overwrites the location of the output directory. Provide the full path.

=head2 lastz_options (def: --notransition --step=30 --seed=match12 --exact=50 --matchcount=1000 --masking=3)

Change the options for lastz. With these default options, lastz run swiftly but will only provide
alignments on regions of very high similarity

Refer to lastz help page for more information.

=head2 axt_chain_options (def: -linearGap=loose)

Change the options for axtChain. "-linearGap=loose" is typically used for chicken/human alignments.
In this application, you probably want to use these to link as many raw lastz alignments as
possible, but you may want to experiment with "-linearGap=medium" or by changing any of the other
options.

Refer to axtChain help page for more information.

=head2 length_threshold (def: 10K)

Change the length threshold for the minimum length of chain to be considered by chain2rdotplot.pl.
This is also used as a maximum gap allowed in any chain before splitting it in two.

=head2 min_length (def: length_threshold)

=head2 max_gap (def: length_threshold)

=head1 REQUIREMENTS

=over

=item lastz

Available from the Miller's lab web page: 

=item axtChain

Part of Jim Kent's libraries and tools: git://genome-source.cse.ucsc.edu/kent.git

=item Rscript

Part of the R package:

=back

=head1 EXAMPLE

  ## ==============================================
  ## 1. Organise your input
  ## ==============================================

  # Create the input directory (called "genomes" by default)

  mkdir genomes

  cd genomes

  # Download the files from Ensembl (for this example, we will use
  # only chr1 and chr2 from human and chimp)

  ../src/download_ensembl_genome.pl --species homo_sapiens --release 77 \
    --chr 1 --chr 2

  ../src/download_ensembl_genome.pl --species pan_troglodytes --release 77 \
    --chr 1 --chr 2A --chr 2B

  gunzip */*.fa.gz

  # Get the length of the chromosomes:

  rm -f *.txt;
  ls -d1 * | while read dir; do perl -lne '\
    if (/^>(\S+)/) {
      if ($name) {
        print "$name\t$l"
      }
      $name = $1;
      $l = 0;
      next
    }
    $l += length($_);
    END {print "$name\t$l"
  }' $dir/*.fa >> $dir.txt; done


  # Run faToNib on each file (naming the nib file according to the
  # sequence name, not the file name)

  ls */*.fa |  perl -lne '
    my $dir = $_;
    $dir =~ s/\/[^\/]*$//;
    my $name = qx"head -n 1 $_";
    $name =~ s/^>(\S+).*/$1/;
    chomp($name);
    print qx"faToNib $_ $dir/$name.nib"'

  cd ..

  ## ==============================================
  ## 2. Configure and run the pipeline
  ## ==============================================

  export PERL5LIB=${PERL5LIB}:src

  init_pipeline.pl src/OxfordPlots.pm --pipeline_url sqlite:///oxford_plot \
    --ref_genome homo_sapiens

  beekeeper.pl --url sqlite:///oxford_plot --loop 

  ## ==============================================
  ## 3. All done!
  ## ==============================================

  ls alignments/*.pdf

=cut


package OxfordPlots;

use strict;
use warnings;
use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');


sub default_options {
    my ($self) = @_;
    return {
        
        %{ $self->SUPER::default_options() },               # inherit other stuff from the base class

        'pipeline_name' => 'oxford_plots',                  # name used by the beekeeper to prefix job names on the farm
#         'host'          => 'tgac-db1',
#         'port'          => '3306',
#         'user'          => 'tgac',
#         'password'      => $self->o('password'),

#         'pipeline_db'   => {
#             -driver => 'mysql',
#             -host   => $self->o('host'),
#             -port   => $self->o('port'),
#             -user   => $self->o('user'),
#             -pass   => $self->o('password'),
#             -dbname => $self->o('ENV', 'USER').'_'.$self->o('pipeline_name'),  # example of a linked definition (resolved via saturation)
#         },
        'pipeline_url'   => 'sqlite:///'.$self->o('pipeline_name'),
       
        'capacity'  => 10,                                 # how many commands can be run in parallel
        
        'work_dir'   => $ENV{PWD},
        'input_dir'  => $self->o('work_dir').'/genomes/',
        'output_dir' => $self->o('work_dir').'/alignments/',
        
        'ref_genome' => $self->o('ref_genome'),

        'lastz_options' => '--notransition --step=30 --seed=match12 --exact=50 --matchcount=1000 --masking=3',
        'axt_chain_options' => '-linearGap=loose',

        'length_threshold' => '1K',
         # ignore any collinear block that is shorter than min_length
        'min_length'       => $self->o('length_threshold'),
         # split any block that has a gap longer than max_gap
        'max_gap'          => $self->o('length_threshold'),
    };
}



sub pipeline_analyses {
    my ($self) = @_;
    
    return [

    ################################################################################################
    ## other_genomes
    ################################################################################################
    {   -logic_name => 'other_genomes',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -parameters => {
            'column_names' => [ 'dir1' ],
        },
        -input_ids => [
        { 
            'input_dir'    => $self->o('input_dir'),
            'inputcmd' => 'ls -1 #input_dir# | grep -v -e .txt -e '.$self->o('ref_genome'), },
        ],
        -flow_into => {
            # will create a fan of jobs
            '2->A' => { 'other_genome_files' => {'input_dir' => '#input_dir#', 'dir1' => '#dir1#'} },   
            'A->2' => { 'make_plot' => {'input_dir' => '#input_dir#', 'dir1' => '#dir1#'} },   
        },
        -meadow_type    => 'LOCAL',

    },
    
    
    ################################################################################################
    ## other_genomes_files
    ################################################################################################
    {   -logic_name => 'other_genome_files',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -parameters => {
            'column_names' => [ 'file1' ],
            'inputcmd' => 'ls -1 #input_dir#/#dir1#/ | grep ".fa$"',
        },
        -flow_into => {
            # will create a fan of jobs
            2 => { 'ref_genome_files' => {'input_dir' => '#input_dir#', 'dir1' => '#dir1#', 
                   'file1' => '#file1#'} },
        },
        -meadow_type    => 'LOCAL',
    },


    ################################################################################################
    ## ref_genome_files
    ################################################################################################
    {   -logic_name => 'ref_genome_files',
        -module     => 'Bio::EnsEMBL::Hive::RunnableDB::JobFactory',
        -parameters => {
            'column_names' => [ 'file2' ],
            'dir2'     => $self->o('ref_genome'),
            'inputcmd' => 'ls -1 '.$self->o('input_dir').'/#dir2#/ | grep ".fa$"',
        },
        -flow_into => {
            # will create a fan of jobs
            2 => { 'lastz' => {'input_dir' => '#input_dir#', 'dir1' =>'#dir1#',
                   'file1' => '#file1#', 'dir2' => '#dir2#', 'file2' => '#file2#'} },
        },
        -meadow_type    => 'LOCAL',
    },
    

    ################################################################################################
    ## lastz
    ################################################################################################
    {   -logic_name     => 'lastz',
        -module         => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters     => {
#            'options'       => '--notransition --step=100 --seed=match12 --exact=250 --matchcount=10000 --masking=3 --chain --format=axt',
#            'options'       => '--notransition --step=30 --seed=match12 --exact=50 --matchcount=1000 --masking=3',

            'options'       => $self->o('lastz_options'),
            'output_dir'    => $self->o('output_dir'),
            'cmd'           => 'mkdir -p #output_dir#/#dir1#; lastz #input_dir#/#dir1#/#file1# #input_dir#/#dir2#/#file2# #options# --format=axt --output=#output_dir#/#dir1#/#file1#.vs.#file2#.axt --rdotplot=#output_dir#/#dir1#/#file1#.vs.#file2#.rdotplot',
        },
        -flow_into => { 1 => "axt_chain" },
    },
    

    ################################################################################################
    ## axt_chain
    ################################################################################################
    {   -logic_name     => 'axt_chain',
        -module         => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters     => {
            'options'       => $self->o('axt_chain_options'),
            'output_dir'    => $self->o('output_dir'),
            'cmd'           => 'axtChain #options# #output_dir#/#dir1#/#file1#.vs.#file2#.axt #input_dir#/#dir1#/ #input_dir#/#dir2# #output_dir#/#dir1#/#file1#.vs.#file2#.chain',
        },
        -flow_into => { 1 => "chain_to_rdotplot" },
    },
    

    ################################################################################################
    ## chain_to_rdoptplot
    ################################################################################################
    {   -logic_name     => 'chain_to_rdotplot',
        -module         => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters     => {
            'output_dir'    => $self->o('output_dir'),
            'min_length'    => $self->o('min_length'),
            'max_gap'    => $self->o('max_gap'),
            'cmd'           => $self->o('work_dir').'/src/chain2rdotplot.pl --chain #output_dir#/#dir1#/#file1#.vs.#file2#.chain --min_length #min_length# --max_gap #max_gap# > #output_dir#/#dir1#/#file1#.vs.#file2#.chain.rdotplot',
        },
    },
    
    
    ################################################################################################
    ## make_plot
    ################################################################################################
    {   -logic_name     => 'make_plot',
        -module         => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
        -parameters     => {
            'output_dir'    => $self->o('output_dir'),
            'cmd'           => 'cd #output_dir#; Rscript '.$self->o('work_dir').'/src/oxford_plot.R #output_dir#/#dir1# #input_dir#; cd ..',
        },
    },
    
    ];
}


1;

