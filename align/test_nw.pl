#!/usr/bin/perl
use strict;
use warnings;
use NW;

tests();

sub tests {
    test_align();
}

sub test_align {
    my $seq1 = "EPSTNL";
    my $seq2 = "PSATL";
    print "\n-------[ test_align ]-------\nSeqs to align:\n1) $seq1\n2) $seq2\n\n";
    my %align_results = NW::align({seq1 => $seq1, seq2 => $seq2, gap => -1});
    my @alignments = @{$align_results{alignments}};
    my $score = $align_results{score};
    print "Results:\n\n";
    print "Score of the alignment: $score\n\n";
    my $i = 1;
    foreach my $alignment (@alignments) {
        print "Alignment number " , $i++, ":\n";
        print "@{$alignment->[0]}", "\n";
        print "@{$alignment->[1]}", "\n";
        print "\n";
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
