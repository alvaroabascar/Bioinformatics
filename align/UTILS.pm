#!/usr/bin/perl
use strict;
use warnings;

package NW;

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
    my @seqs;
    my $file_content;

    # open file and dump content into $file_content
    open (my $fhandle, "<", $filename);
    $file_content = join "", readline $fhandle;
    close($fhandle);
    
    my %seqs = $file_content =~ m/>($seq_names)\n([ATCG\s]+?)\n/g;
    return %seqs;
}

sub load_fasta_seqs {
    return load_fasta_seqs_by_name(shift, "\w\W")
}
