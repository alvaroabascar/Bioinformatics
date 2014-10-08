#!/usr/bin/perl
use strict;
use warnings;

package NW;

#
# Needleman-Wunsch algorithm
#

sub align {
    # align two sequences, return the alignment as an array of two strings
    my %args = %{shift @_};
    my $seq1 = $args{"seq1"};
    my $seq2 = $args{"seq2"};
    my $gap_penalty = $args{"gap_penalty"};
    my $score_matrix = build_score_matrix($seq1, $seq2);
    #DEBUG#
    print "score matrix:\n";
    print_matrix($score_matrix);
    #END_DEBUG
    my $i, my $j;
    my @scores;
    my $arrow_matrix = [];

    # fill column 0 and row 0
    $score_matrix->[0]->[0] = 0;
    $arrow_matrix->[0]->[0] = [];
    for ($i = 1; $i <= length $seq1; $i++) {
        $score_matrix->[$i][0] = $score_matrix->[$i-1][0] + $gap_penalty;
        $arrow_matrix->[$i]->[0] = ["left"];
    }
    for ($j = 1; $j <= length $seq2; $j++) {
        $score_matrix->[0][$j] = $score_matrix->[0][$j-1] + $gap_penalty;
        $arrow_matrix->[0]->[$j] = ["up"];
    }
    # fill 
    for ($i = 1; $i <= length $seq1; $i++) {
        for ($j = 1; $j <= length $seq2; $j++) {
            @scores= ($score_matrix->[$i-1]->[$j-1] + $score_matrix->[$i]->[$j],
                      $score_matrix->[$i-1]->[$j] + $gap_penalty,
                      $score_matrix->[$i]->[$j-1] + $gap_penalty);
            $score_matrix->[$i]->[$j] = max(@scores);
            my @arrows = map { if ($_ == 0) { "diagonal" }
                               elsif ($_ == 1) { "left" }
                               elsif ($_ == 2) { "up" }} max_index(@scores);
            $arrow_matrix->[$i]->[$j] = \@arrows;
        }
    }
    #DEBUG#
    print "\naccumulated score matrix\n";
    print_matrix($score_matrix);
    #END_DEBUG#

    my @seq1_orig = split "", $seq1;
    my @seq2_orig = split "", $seq2;
    return find_alignments({"seq1" => \@seq1_orig,
                                  "seq2" => \@seq2_orig,
                                  "arrow_matrix" => $arrow_matrix});
}

sub find_alignments {
    my %args = %{shift @_};
    my @seq1 = @{$args{"seq1"}};
    my @seq2 = @{$args{"seq2"}};
    my $arrow_matrix = $args{"arrow_matrix"};
    my $i = scalar @seq1;
    my $j = scalar @seq2;
    if ($i == 0 && $j == 0) {
        return (\@seq1, \@seq2);
    }

    # current elements
    my @seq1_aligned = $i > 0 ? ($seq1[$i-1]) : ();
    my @seq2_aligned = $j > 0 ? ($seq2[$j-1]) : ();

    my @result;
    my $i_new;
    my $j_new;
    my @arrows = @{$arrow_matrix->[$i]->[$j]};
    if ("diagonal" ~~ @arrows) {
        $i_new = $i-1;
        $j_new = $j-1;
    } elsif ("left" ~~ @arrows) {
        $i_new = $i-1;
        $j_new = $j;
        shift @seq1_aligned;
        unshift @seq2_aligned, "_";
    } elsif ("up" ~~ @arrows) {
        $i_new = $i;
        $j_new = $j-1;
        shift @seq1_aligned;
        unshift @seq1_aligned, "_";
    }
    my @seq1_new = splice(@seq1, 0, $i_new);
    my @seq2_new = splice(@seq2, 0, $j_new);
    @result = find_alignments({"seq1" => \@seq1_new,
                               "seq2" => \@seq2_new,
                               "arrow_matrix" => $arrow_matrix});
    unshift @seq1_aligned, @{$result[0]};
    unshift @seq2_aligned, @{$result[1]};
    return (\@seq1_aligned, \@seq2_aligned);
}


sub build_score_matrix {
    # build the score matrix, and return a reference to it
    my @seq1 = split "", shift;
    my @seq2 = split "", shift;
    my @score_matrix;

    # first row is made of zeros (due to implementation of align)
    my @first_row = (0) x (scalar @seq2 + 1);
    push @score_matrix, \@first_row;
    foreach my $elem1 (@seq1) {
        # first column is made of zeros
        # put "my @row" here to create a DIFFERENT array for EACH ROW
        my @row = (0);
        foreach my $elem2 (@seq2) {
            push @row, score($elem1, $elem2);
        }
        push @score_matrix, \@row;
    }
    return \@score_matrix;
}

sub score {
    my $elem1 = shift;
    my $elem2 = shift;
    return $elem1 eq $elem2 ? 2 : 0;
}

sub print_matrix {
    my $matrix_ref = shift;
    foreach my $row (@$matrix_ref) {
        foreach my $elem (@$row) {
            print "$elem ";
        }
        print "\n";
    }
}

sub max_index {
    my @arr = @_;
    my $max = shift;
    my @indeces;
    # first find the largest value
    foreach (@_) {
        if ($_ > $max) {
            $max = $_;
        }
    }
    my $i = 0;
    foreach (@arr) {
        if ($_ == $max) {
            push @indeces, $i;
        }
        $i++;
    }
    return @indeces;
}

sub max {
    my $max = shift;
    foreach (@_) {
        if ($_ > $max) {
            $max = $_;
        }
    }
    return $max;
}

1;
