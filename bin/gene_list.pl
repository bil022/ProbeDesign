#!/usr/bin/perl

($genes, $transcript)=@ARGV;
open GENE, "<$genes" or die "gene_lst?";
while (<GENE>) { chomp();
  ($id)=split();
  $G{$id}++;
}

open TRANS, "<$transcript" or die "transcriptW?";
while (<TRANS>) {
  if (/^>/) {
    print STDERR "$gid\t$LEN{$gid}\n" if $gid;
    die unless /gene=(\S+)/;
    ($gid)=($1,$2); $LEN{$gid}=0;
    print if exists $G{$gid}; 
  } elsif (/^([ACTGNactg]+)$/) {
    $LEN{$gid}+=length $1; 
    print uc($_) if exists $G{$gid}; 
  } else { die "$_ not found"; }
}
print STDERR "$gid\t$LEN{$gid}\n" if $gid;

#foreach $gid (keys %LEN) {
#  next unless exists $G{$gid};
#  print STDERR "$gid\t$LEN{$gid}\n";
#}
