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
    my @seq1 = split "", $args{"seq1"};
    my @seq2 = split "", $args{"seq2"};
    my $gap = $args{"gap"};
    my @matrix = build_matrix({seq1 => \@seq1, seq2 => \@seq2,
                               gap => $gap});
    print_matrix(\@matrix);
    return find_alignments({seq1 => \@seq1, seq2 => \@seq2,
                            matrix => \@matrix});
}

sub build_matrix {
    # build a matrix of scores and a matrix of "arrows", return
    # a vector with a reference to each of them.

    my %args = %{shift @_};
    my @seq1 = @{$args{seq1}};
    my @seq2 = @{$args{seq2}};
    my $gap = $args{gap};
    my @matrix;
    # fill upper left element
    $matrix[0][0]{score} = 0;

    my $i;
    my $j;
    # fill first row
    for ($i = 1; $i <= scalar @seq1; $i++) {
        $matrix[$i][0]{score} = $gap * $i;
        $matrix[$i][0]{left} = 1;
    }
    
    # fill first column
    for ($j = 1; $j <= scalar @seq2; $j++) {
        $matrix[0][$j]{score} = $gap * $j;
        $matrix[0][$j]{up} = 1;
    }

    # fill the rest of the matrix
    my $score;
    my @scores;
    my $max_score;
    my @max_scores_indx;
    for ($i = 1; $i <= scalar @seq1; $i++) {
        for ($j = 1; $j <= scalar @seq2; $j++) {
            # score of the pair seq1[i], seq2[j]
            $score = score($seq1[$i-1], $seq2[$j-1]);

            # the three possible scores, corresponding to a diagonal,
            # left and up pointer (arrow)
            @scores = ($score + $matrix[$i-1][$j-1]{score}, # diag
                                $matrix[$i-1][$j]{score} + $gap, # left
                                $matrix[$i][$j-1]{score} + $gap); # up
            $max_score = max(@scores);
            $matrix[$i][$j]{score} = $max_score;

            # index (or indexes, if tie) of the maximal score
            @max_scores_indx = max_index(@scores, $max_score);

            # index 0 = diagonal, index 1 = left, index 2 = up
            my @arrows =  map { ($_ == 0) ? (diagonal => 1) :
                                ($_ == 1) ? (left => 1) :
                                (up => 1) } @max_scores_indx;
                                     
            foreach my $arrow (@arrows) {
                $matrix[$i][$j]{$arrow} = 1;
            }
        }
    }
    return @matrix;
}

sub find_alignments {
    my %args = %{shift @_};
    my @matrix = @{$args{"matrix"}};
    my @seq1 = @{$args{"seq1"}};
    my @seq2 = @{$args{"seq2"}};
    my $i = defined $seq1[0] ? scalar @seq1 : 0;
    my $j = defined $seq2[0] ? scalar @seq2 : 0;
    print "find_alignments $i $j:\n", @seq1, "\n", @seq2, "\n";

    if ($i <= 1 && $j <= 1) {
        print "returning undef\n";
        return ();
    }
    # we align backwards, starting with the last elements of the sequence
    # and then moving towards the begining
    my @seq1_aligned = ($seq1[$i-1]) ;
    my @seq2_aligned = ($seq2[$j-1]);

    my @all_alignments;
    my @seq1_new;
    my @seq2_new;
    my @results;
    if (defined $matrix[$i][$j]->{diagonal}) {
        print "diagonal\n";
        @seq1_new = @seq1[0..$i-2];
        @seq2_new = @seq2[0..$j-2];
        @results = find_alignments({seq1 => \@seq1_new, seq2 => \@seq2_new,
                                    matrix => \@matrix});
        print "results length: ", scalar @results, "\n";
        print "entering loop\n";
        my @alignment;
        if (not @results) {
            push @alignment, [\@seq1_aligned, \@seq2_aligned];
        } else {
            foreach my $result (@results) {
                print "in-loop\n";
                @alignment = @$result;
                print "result: @alignment: $alignment[0] $alignment[1]\n";
                push @{$alignment[0]}, @seq1_aligned;
                push @{$alignment[1]}, @seq2_aligned;
            }
        }
        push @results, @alignment;
    }

    if (defined $matrix[$i][$j]->{up}) {
        print "up\n";
        unshift @seq1_aligned, "_";
        @seq1_new = @seq1[0..$i];
        @seq2_new = @seq2[0..$j-1];
        @results = find_alignments({seq1 => \@seq1_new, seq2 => \@seq2_new,
                                    matrix => \@matrix});
        foreach my $result (@results) {
            if (defined $result) {
                my @alignment = @$result;
                push @{$alignment[0]}, @seq1_aligned;
                push @{$alignment[1]}, @seq2_aligned;
                push @results, @alignment;
            } else {
                push @results, [\@seq1_aligned, \@seq2_aligned];
            }
        }
    }
    if (defined $matrix[$i][$j]->{left}) {
        print "left\n";
        unshift @seq2_aligned, "_";
        @seq1_new = @seq1[0..$i-1];
        @seq2_new = @seq2[0..$j];
        @results = find_alignments({seq1 => \@seq1_new, seq2 => \@seq2_new,
                                    matrix => \@matrix});
        foreach my $result (@results) {
            if (defined $result) {
                my @alignment = @$result;
                push @{$alignment[0]}, @seq1_aligned;
                push @{$alignment[1]}, @seq2_aligned;
                push @results, \@alignment;
            } else {
                push @results, [\@seq1_aligned, \@seq2_aligned];
            }
        }
    }
    print "returning: @results";
    return @results;
}

sub pairs {
    my @arr1 = shift;
    my @arr2 = shift;
    my $len1 = scalar @arr1;
    my $len2 = scalar @arr2;
    my $min = ($len1 < $len2) ? $len1 : $len2;
    my @pairs;
    for (my $i = 0; $i < $min; $i++) {
        push @pairs, [$arr1[$i], $arr2[$i]];
    }
    return \@pairs;
}
sub score {
    my $elem1 = shift;
    my $elem2 = shift;
    return $elem1 eq $elem2 ? 2 : 0;
}

# just for debug..
sub print_matrix {
    my @matrix = @{shift @_};
    my $pad = 5;
    foreach my $row (@matrix) {
        print "|";
        foreach my $elem (@$row) {
            printf "%${pad}d",  $elem->{score}, " ";
        }
        print " " x ($pad - 1), "|\n";
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
