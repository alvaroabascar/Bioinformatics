#!/usr/bin/env perl

use lib './lib';
use strict;
use warnings;
use PWM;

main();

######## STUFF ########

sub main {
  # parse arguments
  my ($inputfile, $n_seqs, $outstream, $log) = parse_args();
  # create PWM object to load and store the PWM, and produce the sequences
  my $pwm = PWM::new('PWM', $inputfile);
  my $seqs = $pwm->n_highest_score_seqs($n_seqs, $log);
  # print seqs
  print_seqs($outstream, $seqs);
}

sub print_seqs {
  my $stream = shift;
  my $sequences = shift;
  foreach my $seq (@$sequences) {
    print $stream "$seq->{seq}\t\t$seq->{score}\n";
  }
}

sub parse_args {

  my ($n_seqs, $inputfile, $outstream, $outputfile, $log);
  my (%errors, @available_opts, $usage, $help);

  $usage = "Usage: main.pl [OPTIONS]... TRANSFAC_FILE\n" .
           "Try './main.pl --help' for more information.\n";

  $help = <<'HELP';
Usage: main.pl [OPTIONS]... TRANSFAC_FILE
Produce the topmost scoring sequences from a position weight matrix.

  -n          number of sequences to be produced (by default, 10)
  -o          output file (by default, STDOUT)
  --log       use log-likelihoods as score (by default, absolute
              frequencies are used)
HELP

  @available_opts = ('-n', '-o', '--score');

  if (not @ARGV) {
    die("Error: missing file operand.\n$usage\n");
  }

  # if asking for help...
  if ($ARGV[0] eq '--help') {
    print($help);
    exit(0);
  }
  # last argument must be the input file
  $inputfile = pop @ARGV;

  # log-likelihood (log = 1) or absolute scores (log = 0)?
  $log = 0;
  my $pos;
  foreach my $i (0..$#ARGV) {
    if ($ARGV[$i] eq '--log') {
      $log = 1;
      $pos = $i;
      last;
    }
  }
  if ($log) {
    splice @ARGV, $pos, 1;
  }
  # the other arguments are "-option argument" => even number of tokens
  if (scalar @ARGV % 2 != 0) {
    die("Error:\n$usage");
  }
  # turn args into a dictionary option => argument
  my %args = @ARGV;

  # error if some option is unknown
  foreach my $key (keys %args) {
    if (not grep { $_ eq $key } @available_opts) {
      die("Unknown option: $key\n$usage");
    }
  }

  # number of sequences is 10 by default
  if ($args{'-n'}) {
    $n_seqs = $args{'-n'};
  }
  else {
    $n_seqs = 10;
  }

  # outputfile is optional. If not provided, use STDOUT as output stream
  $outputfile = $args{'-o'};
  if ($outputfile) {
    open($outstream, ">$outputfile") || die $!;
  }
  else {
    $outstream = *STDOUT;
  }

  # an informative message
  if ($args{'-o'}) {
    print "Output saved to \"$args{'-o'}\"\n";
  }
  return ($inputfile, $n_seqs, $outstream, $log);
}

