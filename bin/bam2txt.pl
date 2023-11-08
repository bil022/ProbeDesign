#!/usr/bin/perl

while (<>) {
  ($id, $flag, $chr, $s, $cigar, $seq)=split();
  die $seq unless $seq=~/^[ACTG]+$/;
  $revcomp = reverse $seq;
  $revcomp =~ tr/ATGCatgc/TACGtacg/;
  $rev=$flag&16;
  $key="$chr:$s:$flag:$cigar:$seq";
  $seq=$revcomp unless $rev; # reverse if plus strand
  print "$id\t$seq\t$key\n";
}

__END__
ECE1.19.0.28	16	chr1	21217369	30M	AGCCAGCCCTATGGCTCCCCAGGCCTAAGG
