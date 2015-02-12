#!/usr/bin/env perl
use strict;
use warnings;
use Seq;
use Data::Dumper;

package PWM;

sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;
  $self->load_from_file(shift) if (@_);
  $self->{log} = 0;
  return $self;
}

sub load_from_file {
  # Load a PWM matrix in TRANSFAC format from a file, return a dictionary:
  #  { id => motif name
  #    bf => species name
  #    length => length of the motif
  #    alphabet => ordered alphabet (might be nucleotides, aminoacids...)
  #    matrix => [
  #               [matrixA, matrixT, counC, matrixG], # row1
  #               ... #..till row length-1
  #               ]
  #  }
  #
  # (This function works for any type of sequence)

  my $self = shift;
  my $filename = shift;
  open (my $fh, "<$filename") or die $!;
  my @cols;
  my @matrix;
  $self->{length} = 0;
  while (<$fh>) {
    chomp $_;
    # Split line into columns
    @cols = split " ", $_;

    # ID:
    if ($cols[0] eq "ID") {
      $self->{ID} = $cols[1];

    # No ID yet? keep searching
    } elsif (not defined $self->{ID}) {
      next;

    # Counts:
    } elsif ($cols[0] =~ m/^[0-9]+/) {
      my %row;
      foreach my $i (1..$#cols - 1) {
        $row{$self->{alphabet}->[$i-1]} = $cols[$i];
      }
      $row{order} = $self->{alphabet};
      push @matrix, \%row;
      $self->{length}++;

    # Species name:
    } elsif ($cols[0] eq "BF") {
      $self->{BF} = join " ", @cols[1..$#cols];

    # Alphabet
    } elsif ($cols[0] eq "P0") {
      $self->{alphabet} = [@cols[1..$#cols]];

    # End (//)
    } elsif ($cols[0] eq "//") {
      last;
    }
  }
  close($fh);
  $self->{matrix} = \@matrix;
  $self->order_rows_by_score();
}

sub order_rows_by_score {
  # For each row in the given matix, reorder the columns in decr. order
  # and set as self->{matrix} this new matrix:
  # what changes of order is the vector $row->{order}

  my $self = shift;
  my $matrix = $self->{matrix};
  for my $row_ref (@$matrix) {
    # in $self->matrix each row is a hash { letters => scores}
    # now each row is an array of hashes { score => X, letter => Y}
    # (an array can be ordered, but not a hash)
    # Add row to matrix, ordered from higher to lower score
    $row_ref->{order} = [sort { $row_ref->{$b} <=> $row_ref->{$a} }
                              @{$row_ref->{order}}];
  }
}

sub get_highest_score_seq {
  # Given the Position Weight Matrix, return the
  # sequence with the highest score (any of them in case of tie)
  my $self = shift;
  # elem 0 corresponds to the highest scoring letter of each position
  my $highest_seq = [map { 0 } @{$self->{matrix}}];
  my $seq = Seq::new('Seq', $highest_seq, 4, $self);
  return $seq;
}

sub change_size {
  # Given a row and two indexes, computes the difference of scores among
  # the score at matrix[row][index1] and matrix[row][index2]
  my $self = shift;
  my $args = shift;
  my $row = $args->{row};
  my $index1 = $args->{index1};
  my $index2 = $args->{index2};
  my $change_size = $self->get_score_by_index({row   => $row,
                                               index => $index1 })
                  - $self->get_score_by_index({row   => $row,
                                               index => $index2 });
  return $change_size;
}

sub get_elem_by_index {
  # Retrieve the element at row $row, which is the $n(th) biggest
  my $self = shift;
  my $args = shift;
  my $letter = $self->{matrix}[$args->{row}]{order}[$args->{n}];
  return $letter;
}

sub get_score_by_index {
  # Retrieve the score at row $row, which is the $n(th) biggest
  my $self = shift;
  my $args = shift;
  my $matrix = $self->{matrix};
  my $letter = $matrix->[$args->{row}]->{order}->[$args->{index}];
  my $score = $matrix->[$args->{row}]->{$letter};
  return $score;
}

sub get_index_by_elem {
  # Given a row and an element, return it's position (ie. it's index in the
  # ordered row)
  my $self = shift;
  my $args = shift;
  my $elems = $self->{matrix}->[$args->{row}]->{order};
  my ($index) = grep { $elems->[$_] eq $args->{elem} } 0..$#{$elems};
  return $index;
}

sub get_score_by_elem {
  # Retrieve score at row args->row, for the element args->elem
  my $self = shift;
  my $args = shift;
  my $row = $args->{row};
  my $elem = $args->{elem};
  my $matrix = $self->{matrix};
  my $score = $matrix->[$row]->{$elem};
  return $score;
}

sub freqs_sum_in_row {
  # Return the sum of frequencies in a given row
  my $self = shift;
  my $args = shift;
  my $row = $args->{row};
  my $matrix = $self->{matrix};
  my $sum = 0;
  foreach my $nucl (@{$self->{alphabet}}) {
    $sum += $matrix->[$row]->{$nucl};
  }
  return $sum;
}

sub score {
  # Given a sequence, calculate the score
  my $self = shift;
  my @seq = @{shift @_};
  my $log = $self->{log};
  return undef unless ($#seq eq $#{$self->{matrix}});
  my $score = 0;
  my ($itemfreq, $totalfreq);
  # compute score as the sum of absolute frequencies
  if (not $log) {
    foreach my $i (0..$#seq) {
      $score += $self->get_score_by_index({ row => $i, index => $seq[$i] });
    }
  }
  # compute score as the sum of log likelihoods
  else {
    foreach my $i (0..$#seq) {
      $itemfreq = $self->get_score_by_index({ row => $i, index => $seq[$i] });
      $totalfreq = $self->freqs_sum_in_row({row => $i});
      $itemfreq = $itemfreq eq 0 ? 0.01 : $itemfreq;
      $score += log($itemfreq / $totalfreq);
    }
  }
  return $score;
}

sub n_highest_score_seqs {
  # Given the PWM, return the "n" seqs with the highest score, ordered from
  # higher to lower score
  my $self = shift;
  my $n = shift;
  # 1 -> return scores as log likelihoods. 0 -> return them as sum of freqs
  my $log = shift;
  $self->{log} = $log;
  my @matrix = @{$self->{matrix}};

  # A container with the produced sequences, which are expected to be
  # modified to produce new sequences with a high score
  my @seqs_bag = ();

  # STEP 1 of the algorithm: obtain the highest score sequence

  my $highest_seq = $self->get_highest_score_seq();
  push @seqs_bag, $highest_seq;

  my ($best_new_score, @smallest_changes, $change_size, $new_score);
  my $i = 0;
  my $seqs_so_far = 1;
  while ($seqs_so_far < $n) {
    # Here we store the best new score found so far
    $best_new_score = 0;

    # Here we store the smallest changes, with the format (seq => , pos =>)
    @smallest_changes = ();

    # for each high scoring sequence already found...
    foreach my $seq (0..$#seqs_bag) {
      my $currseq = $seqs_bag[$seq];

      # if this sequence has already been modified in all positions,
      # skip it to avoid producing a repeated sequence
      if ($currseq->all_unmodifiable()) {
        next;
      }

      # for each position in the sequence...
      foreach my $pos (0..$currseq->length()) {

        # if this position had already been modified (and thus, is
        # unmodifiable), or if the element at that position is the lowest
        # possible, then jump to next position
        if ((not $currseq->modifiable($pos)) or
            ($currseq->element($pos) + 1 >= scalar @{$self->{alphabet}})) {
          next;
        }
        $change_size = $self->change_size({row => $pos,
                                           index1 => $currseq->element($pos),
                                           index2 => $currseq->element($pos)+1
                                          });
        $new_score = $seqs_bag[$seq]->{score} - $change_size;

        if ($new_score == $best_new_score) {
          push @smallest_changes, {seq => $seq, pos => $pos};
        }
        elsif (abs $new_score > abs $best_new_score) {
          $best_new_score = $new_score;
          @smallest_changes = ({seq => $seq, pos => $pos});
        }
      }
    }

    # haven't found any good change?
    if (not @smallest_changes) {
      return $self->return_seqs({ seqs => \@seqs_bag});
    }
    # by now we have found the smallest change(s)
    foreach my $change (@smallest_changes) {
      if ($seqs_so_far >= $n) {
        last;
      }
      my ($seq, $pos) = ($change->{seq}, $change->{pos});

      my $foundseq = $seqs_bag[$seq];
      # set that position of the sequence as modified, in order to avoid
      # modifying it again
      $foundseq->set_unmodifiable($pos);

      # Create a new sequence:
      # 1. clone the selected sequence
      my $newseq = $foundseq->clone();

      # 2. apply the change to $pos (increment by one)
      $newseq->change($pos, $newseq->element($pos) + 1);

      # this avoids repetition, but multiplies executing time by more than 2...
      if ($self->already_in($newseq, \@seqs_bag)) {
        next;
      }
      # if the new index is the last possible, don't try to modify it
      # later
      if ($newseq->element($pos) == $#{$self->{alphabet}}) {
        $newseq->set_unmodifiable($pos);
      }

      # 3. By now it has the score of its parent. Update it!
      $newseq->set_score($self->score($newseq->{seq}));

      # 4. Add it to seqsbag
      push @seqs_bag, $newseq;

      $seqs_so_far++;
    }
  }
  return $self->return_seqs({ seqs => \@seqs_bag });
}

#
# Conversion of a seq of indices to real seq of letters, according
# to the Position Weight Matrix
#
sub seq_to_realseq {
  my $self = shift;
  my $seq = shift;
  my $matrix = $self->{matrix};
  my @realseq = ();
  my $i = 0;
  return [map { $matrix->[$i++]->{order}->[$_] } @{$seq->sequence()}];
}

sub return_seqs {
  my $self = shift;
  my $args = shift;
  my @seqs = @{$args->{seqs}};
  return [map { { seq => join("", @{$self->seq_to_realseq($_)}),
                         score => $_->{score}} } @seqs];
}

sub print_matrix {
  my $self = shift;
  my $matrix = $self->{matrix};
  foreach my $ref (@$matrix) {
    foreach my $elem (@{$ref->{order}}) {
      print STDERR $ref->{$elem}, " ";
    }
    print STDERR "\n";
  }
}

sub already_in {
  my $self = shift;
  my $myseq = shift;
  my $arr_seqs = shift;
  my $j1 = join "", @{$myseq->{seq}};
  foreach my $seq (@$arr_seqs) {
    my $j2 = join "", @{$seq->{seq}};
    if ($j1 eq $j2) {
      return 1;
    }
  }
  return 0;
}

return(1);
