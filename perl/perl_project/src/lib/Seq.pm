#!/bin/usr/env perl
use strict;
use warnings;

package Seq;
# The class Seq represents a sequence, with some additional methods

sub new {
  my $class = shift;
  my $seq = shift;
  my $alphabet_length = shift;
  my $weight_matrix = @_ ? shift : undef;
  # ref to an array of indexes
  # all pos modifiable
  my $self = { seq => $seq,
               modifiable => [map { $_ < $alphabet_length-2 ? 1 : 0 } @$seq],
               score => $weight_matrix ? $weight_matrix->score($seq) : undef,
               alphabet_length => $alphabet_length
             };
  return bless $self, $class;
}

# Retrieve the full sequence
sub sequence {
  my $self = shift;
  return $self->{seq};
}

# Retrieve whether the position $pos can be modified
sub modifiable {
  my $self = shift;
  my $pos = shift;
  return $self->{modifiable}->[$pos];
}

# Retrieve the element at pos $pos
sub element {
  my $self = shift;
  my $pos = shift;
  return $self->{seq}->[$pos];
}
#
# Return a copy of this sequence object
#
    # THIS FUNCTION IS KEY TO PROGRAM PERFORMANCE #
sub clone {
  my $self = shift;
  return new('Seq', [@{$self->{seq}}], $self->{alphabet_length});

  #### 
  # don't do this!!!! you can't use the same reference... 
  # return new('Seq', $self->{seq}, $self->{alphabet_length});
  ####
}

# Return the length of the sequence
sub length {
  my $self = shift;
  return $#{$self->{seq}};
}

# set position $pos as unmodifiable
sub set_unmodifiable {
  my $self = shift;
  my $pos = shift;
  $self->{modifiable}->[$pos] = 0;
}

# Change the element at pos $pos by $newelem
sub change {
  my $self = shift;
  my $pos = shift;
  my $newelem = shift;
  $self->{seq}->[$pos] = $newelem;
}

# Return 1 if all positions are unmodifiable, else 0
sub all_unmodifiable {
  my $self = shift;
  return (not scalar grep { $_ eq 1 } @{$self->{modifiable}});
}

sub set_score {
  my $self = shift;
  my $score = shift;
  $self->{score} = $score;
}

return(1);
