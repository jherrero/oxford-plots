#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Long;

my $help;
my $verbose = 0;
my $chain_file;
my $min_length = 10000000;
my $max_gap = 10000000;

my $desc =qq{SCRIPT
    
    chain2rdotplot.pl
    
DESCRIPTION:
    
    This script transforms an axtChain file into an rdotplot file.
    
SYNOPSIS:

    chain2rdotplot.pl [options] file.chain > file.rdotplot
    
OPTIONS:
    
    -help
        Shows this help
    
    -distance_threshold distance
        This defines the maximum gap allowed within a chain and
        the minimum length of a chain
        Def:

    -chain filename
        Location of the chain file (output of axtChain). The first
        unnamed argument is used if not defined.

PLOTTING WITH R:
    
    dots <- read.table("file.rdotplot");
    lines(dots);

AUTHOR:

    Javier Herrero (Javier.Herrero\@tgac.ac.uk)

};

GetOptions(
    "help" => \$help,
    "verbose" => \$verbose,
    "chain=s" => \$chain_file,
    "min_length=s" => \$min_length,
    "max_gap=s" => \$max_gap,
);

if ($help) {
    print $desc;
    exit(0);
}

if ($min_length =~ /^(\d+)M$/) {
    $min_length = $1 * 1000000;
} elsif ($min_length =~ /^(\d+)K$/) {
    $min_length = $1 * 1000;
}

if ($max_gap =~ /^(\d+)M$/) {
    $max_gap = $1 * 1000000;
} elsif ($max_gap =~ /^(\d+)K$/) {
    $max_gap = $1 * 1000;
}

if (!$chain_file and @ARGV) {
    $chain_file = shift(@ARGV);
}

open(CHAIN, $chain_file) or die;

my $this_chain = undef;
my $is_header_printed = 0;

while (<CHAIN>) {
    if (/^chain/) {
        
        print_chain($this_chain) if ($this_chain);

        my ($seq1, $len1, $str1, $start1, $end1, $seq2, $len2, $str2, $start2, $end2) = $_ =~ /chain \d+ (\S+) (\d+) (.) (\d+) (\d+) (\S+) (\d+) (.) (\d+) (\d+)/;
        if (!$is_header_printed) {
            print "$seq1\t$seq2\n";
            $is_header_printed = 1;
        }
        
        if ($str1 eq "-") {
            my $temp = $len1 - $start1;
            $end1 = $len1 - $start1;
            $start1 = $temp;
        }
        if ($str2 eq "-") {
            my $temp = $len2 - $start2;
            $end2 = $len2 - $end2;
            $start2 = $temp;
        }

        $this_chain = undef;
        
        next if (abs($end1 - $start1) < $min_length);
        next if (abs($end2 - $start2) < $min_length);

        $this_chain = [$start1, $start1, ($str1 eq "+")?1:-1, $start2, $start2, ($str2 eq "+")?1:-1];

    } elsif ($this_chain and /^(\d+)\t(\d+)\t(\d+)$/) {
        if ($2 > $max_gap or $3 > $max_gap) {
            print_chain($this_chain);
            $this_chain->[0] = $this_chain->[1] + $this_chain->[2] * ($1 + $2);
            $this_chain->[3] = $this_chain->[4] + $this_chain->[5] * ($1 + $3);
        }
        $this_chain->[1] += $this_chain->[2] * ($1 + $2);
        $this_chain->[4] += $this_chain->[5] * ($1 + $3);
    
    } elsif ($this_chain and /^(\d+)$/) {
        $this_chain->[1] += $this_chain->[2] * $1;
        $this_chain->[4] += $this_chain->[5] * $1;
    }
}
print_chain($this_chain) if ($this_chain);

sub print_chain {
    my ($this_chain) = @_;

    if (abs($this_chain->[1] - $this_chain->[0]) < $min_length) {
        return;
    }
    if (abs($this_chain->[4] - $this_chain->[3]) < $min_length) {
        return;
    }

#    print join("\t", @$this_chain), "\n";
    
    print join("\n",
        $this_chain->[0]."\t".$this_chain->[3],
        $this_chain->[1]."\t".$this_chain->[4],
        "NA\tNA"), "\n";
}
