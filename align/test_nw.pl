#!/usr/bin/perl
use strict;
use warnings;
use NW;

tests();

sub tests {
    test_align();
}

sub test_align {
    #my $seq1 = "JUPSAL";
    #my $seq2 = "EPSTNL";
    my $seq1 = "JUPSAL";
    my $seq2 = "EUPSAL";
    print "\ntest_align. Aligning:\n$seq1\n$seq2\n\nResult:\n\n";
    my @alignments = NW::align({seq1 => $seq1, seq2 => $seq2, gap => -1});
    print "len result:\n", scalar @alignments;
    foreach my $alignment (@alignments) {
        print $alignment->[0], "\n";
        print $alignment->[1], "\n";
    }
}

sub test_build_matrix {
    my $seq1 = "GCATGCU";
    my $seq2 = "GATTACA";
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
