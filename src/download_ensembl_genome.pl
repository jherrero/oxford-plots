#! /usr/bin/env perl
use strict;
use warnings;

use Net::FTP;
use File::Path qw(make_path);
use Getopt::Long;

my $help = 0;
my $release = 73;
my $species;
my @chrs;
my $verbose = 0;

GetOptions(
    "help"  => \$help,
    "species=s" => \$species,
    "release=i" => \$release,
    "chr=s@" => \@chrs,
    "verbose" => \$verbose,
    );

if ($help or !$species) {
    print "USE: download_ensembl_genome.pl --species homo_sapiens --release 73 [options]\n",
          " where options can be --chr 1 --chr X\n";
    exit(0);
}

download_genome($species, $release);

sub download_genome {
    my ($species, $release) = @_;

    make_path("$species");
    chdir("$species");
    my $ftp=Net::FTP->new("ftp.ensembl.org",Debug=>0);
    $ftp->login("anonymous","anonymous");
    $ftp->binary();
    print STDERR "Moving to dir /pub/release-$release/fasta/$species/dna/...\n";
    my $ftp_dir = "/pub/release-$release/fasta/$species/dna/";
    $ftp->cwd($ftp_dir);
    my @files = $ftp->ls();

    foreach my $this_file (grep {/\.dna_sm\.chromosome\.[XY\d].?\.fa.gz/} @files) {
        next if (@chrs and !grep {$this_file =~ /chromosome\.$_\.fa/} @chrs);
        if ($verbose) {
            print STDERR "Downloading $this_file...\n";
        }
        $ftp->get($this_file) or die $ftp->message;
    }
    chdir("..");
}