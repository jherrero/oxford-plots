OxfordPlots
===========

This software is intended to build pairwise dotplots (also known as Oxford plots) for a set of
species vs one particular reference genome. It can be used for just a pair of species or just a
sub-set of chromosomes (instead of using whole genomes).

Usage
=====

This has been developed as an eHive pipeline (see https://github.com/Ensembl/ensembl-hive). Firstly,
you need to download the (soft-masked) sequences, one chromosome per file, one species per
directory. You also need to pre-calculate the length of each chromosome and store this in a text
file.

Please refer to the PerlDoc of src/OxfordPlots.pm for furhter details on the options.

Full Example
============

Organise your input
-------------------

```
# Create the input directory (called "genomes" by default)

mkdir genomes

cd genomes

# Download the files from Ensembl (for this example, we will use only chr1 and chr2 from human and chimp)

../src/download_ensembl_genome.pl --species homo_sapiens --release 77 --chr 1 --chr 2

../src/download_ensembl_genome.pl --species pan_troglodytes --release 77 --chr 1 --chr 2A --chr 2B

gunzip */*.fa.gz

# Get the length of the chromosomes:

rm -f *.txt;
ls -d1 * | while read dir; do perl -lne '
  if (/^>(\S+)/) {
    if ($name) {
      print "$name\t$l"
    };
    $name = $1;
    $l = 0;
    next
  } $l += length($_);
  END {print "$name\t$l"}
  ' $dir/*.fa >> $dir.txt; done

# Run faToNib on each file (naming the nib file according to the sequence name, not the file name)

ls */*.fa |  perl -lne '
 my $dir = $_;
 $dir =~ s/\/[^\/]*$//;
 my $name = qx"head -n 1 $_";
 $name =~ s/^>(\S+).*/$1/;
 chomp($name);
 print qx"faToNib $_ $dir/$name.nib"'

cd ..
```

Configure and run the pipeline
------------------------------

```
export PERL5LIB=${PERL5LIB}:src

init_pipeline.pl src/OxfordPlots.pm --pipeline_url sqlite:///oxford_plot --ref_genome homo_sapiens

beekeeper.pl --url sqlite:///oxford_plot --loop
```

File organisation
=================

(INPUT_DIR after downloading the sequences)
genomes
genomes/homo_sapiens
genomes/homo_sapiens/Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa
genomes/homo_sapiens/Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa
genomes/pan_troglodytes
genomes/pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.1.fa
genomes/pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2A.fa
genomes/pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2B.fa


(INPUT_DIR after calculating the sequence lengths)
genomes
genomes/homo_sapiens
genomes/homo_sapiens/Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa
genomes/homo_sapiens/Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa
genomes/homo_sapiens.txt
genomes/pan_troglodytes
genomes/pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.1.fa
genomes/pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2A.fa
genomes/pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2B.fa
genomes/pan_troglodytes.txt


(INPUT_DIR after running faToNib)
genomes
genomes/homo_sapiens
genomes/homo_sapiens/1.nib
genomes/homo_sapiens/2.nib
genomes/homo_sapiens/Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa
genomes/homo_sapiens/Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa
genomes/homo_sapiens.txt
genomes/pan_troglodytes
genomes/pan_troglodytes/1.nib
genomes/pan_troglodytes/2A.nib
genomes/pan_troglodytes/2B.nib
genomes/pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.1.fa
genomes/pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2A.fa
genomes/pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2B.fa
genomes/pan_troglodytes.txt


(OUTPUT_DIR after running the pipeline)
alignments/
alignments//pan_troglodytes
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.1.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa.axt
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.1.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa.chain
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.1.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa.chain.rdotplot
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.1.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa.rdotplot
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.1.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa.axt
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.1.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa.chain
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.1.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa.chain.rdotplot
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.1.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa.rdotplot
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2A.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa.axt
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2A.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa.chain
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2A.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa.chain.rdotplot
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2A.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa.rdotplot
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2A.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa.axt
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2A.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa.chain
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2A.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa.chain.rdotplot
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2A.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa.rdotplot
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2B.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa.axt
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2B.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa.chain
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2B.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa.chain.rdotplot
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2B.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa.rdotplot
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2B.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa.axt
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2B.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa.chain
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2B.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa.chain.rdotplot
alignments//pan_troglodytes/Pan_troglodytes.CHIMP2.1.4.dna_sm.chromosome.2B.fa.vs.Homo_sapiens.GRCh38.dna_sm.chromosome.2.fa.rdotplot
alignments//Pan_troglodytes.pdf
