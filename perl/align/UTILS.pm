#!/usr/bin/perl
use strict;
use warnings;

package UTILS;

#
# Utilities to work with sequences
#

# We could just load all sequences and then discard the ones we don't need.
# However, matching only the desired sequences is more memory-efficient
# and this might be useful if there are many or very long sequences
sub load_fasta_seqs_by_name {
    # load a the specified fasta sequences from the specified them
    # return them as a hash {name} -> sequence
    my $filename = shift;
    my $seq_names = join "|", @_;
    my $file_content;

    # open file and dump content into $file_content
    open (my $fhandle, "<", $filename);
    $file_content = join "", readline $fhandle;
    close($fhandle);
    
    #my @seqs = $file_content =~ m/>($seq_names)\s+?([ATCG\s]+)\n/g;
    my @seqs = $file_content =~ m/>($seq_names)\s+?([A-Za-z\s]+)\n/g;
    map { $_ =~ s/\s//g } @seqs;
    return @seqs;
}

sub load_fasta_seqs {
    return load_fasta_seqs_by_name(shift, '\w\W');
}

sub pretty_align {
    my @seq1 = @{shift @_};
    my @seq2 = @{shift @_};
    my $cols = join "", `tput cols`;   # "backticks" execute a command and
                                       # return a list of lines as output
    my $i = 0;
    my @compare = map { ($_ eq $seq2[$i++]) ? "*" : "." } @seq1; 
    while (@seq1 || @seq2) {
        print join "", splice(@seq1, 0, $cols-1), "\n";
        print join "", splice(@seq2, 0, $cols-1), "\n";
        print join "", splice(@compare, 0, $cols-1), "\n\n";
    }
}
 
1;   
