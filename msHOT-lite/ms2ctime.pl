#!/usr/bin/perl -w

use strict;
use warnings;
use NHX;
use Getopt::Std;

my %opts = (d=>'ctime', u=>2.5e-8);
die("Usage: ms2ctime.pl [-u $opts{u}] [-d $opts{d}] <in.ms>\n") if (-t STDIN && @ARGV == 0);

my (%data, $nind, $k, $L, @fh, $chr, $N0);
while (<>) {
  if (/^\S*ms\S+\s+(\d+)\s+(\d+)/) { # the command line
	$data{NSAM} = $1; $data{NREP} = $2;
	$data{T} = $1 if (/-t\s+(\S+)/);
	$data{R} = $1, $data{L} = $2 if (/-r\s+(\S+)\s+(\d+)/);
	$data{LH3} = /-l/? 1 : 0;
	$data{N0} = $N0 = $data{T} / (4.0 * $opts{u} * $data{L});
	$chr = 0;
	$nind = int($data{NSAM}/2); # number of diploid individuals
	mkdir($opts{d}) unless (-d $opts{d});
	open($fh[$_], ">$opts{d}/ind-".($_*2-1)."-".($_*2)) for (1 .. $nind);
  } elsif (/^\/\//) {
	$k = $L = 0; ++$chr;
	delete($data{SEG});
  } elsif (/^\[(\d+)\](\(.*\);)/) { # marginal tree
	my $nhx = NHX->new;
	my $seg = \@{$data{SEG}[$k++]};
	$seg->[0] = $1;
	$nhx->parse($2);
	$seg->[$_] = int(get_dist2($nhx, $_*2-1, $_*2) * 2 * $N0 + 0.5) for (1 .. $nind); # 2 == 1/2*4
	$L += $seg->[0];
	if ($L == $data{L}) { # write ctime
	  for my $i (1 .. $nind) {
		my ($last, $begin, $end) = (-1.0, 0, 0);
		for my $p (@{$data{SEG}}) {
		  my $t = $p->[$i];
		  if ($t == $last) {
			$end += $p->[0];
		  } else {
			if ($last >= 0.0) {
			  print {$fh[$i]} "chr$chr\t", $begin+1, "\t$end\t$last\n";
			}
			$last = $t; $begin = $end; $end += $p->[0];
		  }
		} # for $p
		print {$fh[$i]} "chr$chr\t", $begin+1, "\t$end\t$last\n";
	  }
	}
  }
}
close($fh[$_]) for (1 .. $nind);

sub get_dist2 {
  my ($nhx, $p, $q) = @_;
  my $dist = 0;
  $p = $nhx->get_leaf_by_name($p);
  $q = $nhx->get_leaf_by_name($q);
  while ($p != $q) {
	if ($p->{_index} < $q->{_index}) { $dist += $p->{_L}; $p = $p->{_P}; }
	else { $dist += $q->{_L}; $q = $q->{_P}; }
  }
  return $dist;
}
