#!/usr/bin/perl
use strict;
use warnings;
use FindBin;                 # locate this script
use lib "$FindBin::Bin/..";  # use the parent directory
use NW;
use UTILS;

tests();

sub tests {
    test_align_from_file();
}

sub test_align_from_file {
    my $file = "example.fasta";
    my %seqs = UTILS::load_fasta_seqs_by_name($file, "t111", "1pdy");
    test_align(values %seqs);
}

sub test_align {
    my $seq1 = shift;
    my $seq2 = shift;
    print "\n-------[ test_align ]-------\nSeqs to align:\n1) $seq1\n2) $seq2\n\n";
    my %align_results = NW::align({seq1 => $seq1, seq2 => $seq2, gap => -1});
    my @alignments = @{$align_results{alignments}};
    my $score = $align_results{score};
    print "Results:\n\n";
    print "Score of the alignment: $score\n\n";
    my $i = 1;
    foreach my $alignment (@alignments) {
        print "\nAlignment number " , $i++, ":\n";
        UTILS::pretty_align(\@{$alignment->[0]}, \@{$alignment->[1]});
    }
}

sub test_build_matrix {
    my $seq1 = "GAT";
    my $seq2 = "GT";
    say("building a score matrix for:\n$seq1 and\n$seq2\n");
    my @matrix = NW::build_matrix($seq1, $seq2);
    NW::print_matrix(@matrix);
}


sub test_print_matrix {
    my @matrix = [[1, 2, 3],
                  [4, 5, 6],
                  [7, 8, 9]];
    NW::print_matrix(\@matrix);
}

sub say {
    push @_, "\n";
    print @_;
}
