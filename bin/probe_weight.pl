#!/usr/bin/perl

die "Useage: $0 <gene> <gc> <len> <strand>[@ARGV]" unless @ARGV==4;
($gid, $N, $len, $strand)=@ARGV;

sub byProbe {
  @As=split(/\t/, $a); 
  if ($As[0]=~/^(\S+)\.(\d+)\.(\d+)\.(\d+)\.(\d+)$/) {
    ($gA,$gcA,$offA,$qA,$ordA)=($1,$2,$3,$4,$5);
  } elsif ($As[0]=~/^(\S+)\.(\d+)\.(\d+)\.(\d+)$/) {
    ($gA,$gcA,$offA,$ordA)=($1,$2,$3,$4); $qA='1';
  } else { die $_; }

  @Bs=split(/\t/, $b);
  if ($Bs[0]=~/^(\S+)\.(\d+)\.(\d+)\.(\d+)\.(\d+)$/) {
    ($gB,$gcB,$offB,$qB,$ordB)=($1,$2,$3,$4,$5); 
  } elsif ($Bs[0]=~/^(\S+)\.(\d+)\.(\d+)\.(\d+)$/) {
    ($gB,$gcB,$offB,$ordB)=($1,$2,$3,$4); $qB='1';
  } else { die $_; }

  # ordA/B ordered by position, #total=$POS 
  return ord($qB)<=>ord($qA) if ord($qB)!=ord($qA);
  die "$As[0]?$gA!=$gB" unless $gA eq $gB;
  $absA=abs($gcA-$N); $absB=abs($gcB-$N);
  return $absA<=>$absB if $absA!=$absB; 
  return $offA<=>$offB if $offA!=$offB;
  $posA=$As[3]; $posB=$Bs[3];
  die "$a\n$b" if $posA==$posB;
  $strandA=$As[1]&16; $strandB=$Bs[1]&16;
  die "$a\n$b" unless ($strandA == $strandB);
  return $posA<=>$posB if $strandA!=0;
  return $posB<=>$posA; 
}

$POS=1;
while (<STDIN>) { chomp();
  next unless /^$gid\./;
  @data=split(); 
  $flag=$data[1]; $rev=$flag&16;
  if ($rev) { next unless $strand; }
  else { next if $strand; } 
  # print "$flag:$rev:$strand\t$_\n";
  $pos=$data[3];
  next if $pos == $lastPos;
  s/\t/.$POS\t/;
  $CNT{$g}++;
  push(@LNs, $_);
  $lastPos=$pos;
  $POS++;
}

$LN=1;
foreach (sort byProbe @LNs) {
  @data=split();
  # weight w. strand #probes:$CNT{$g} 
  die $_ unless $data[0]=~/^(\S+)\.(\d+)\.(\d+)\.(\d+)$/; ($g,$gc,$numOff,$ord)=($1,$2,$3,$4);
  if ($data[0]=~/^(\S+)\.(\d+)\.(\d+)\.(\d+)$/) {
    ($g,$gc,$numOff,$ord)=($1,$2,$3,$4);
  } elsif ($data[0]=~/^(\S+)\.(\d+)\.(\d+)\.(\d+)\.(.)$/) {
    ($g,$gc,$numOff,$ord,$qual)=($1,$2,$3,$4,$5);
  } else { die $_; }
  # $revStrand=$data[1]&16; 

  s/\t/.$LN\t/;
  print; print "\n";
  $LN++;
}

__END__
LEPR.21.7	0	chr1	65420686	255	30M	*	0	0	CAGGAAGCCGGAAGCAGCCGCGGCCCCAGT	*	NH:i:1	HI:i:1	AS:i:29	nM:i:0
LEPR.19.10	0	chr1	65420708	255	30M	*	0	0	GCCCCAGTTCGGGAGACATGGCGGGCGTTA	*	NH:i:1	HI:i:1	AS:i:29	nM:i:0
LEPR.17.7	0	chr1	65425301	255	30M	*	0	0	AGCTCTCGTGGCATTATCCTTCAGTGGGGC	*	NH:i:1	HI:i:1	AS:i:29	nM:i:0
LEPR.15.4	0	chr1	65425350	255	30M	*	0	0	ATGCTGGGATGTGCCTTAGAGGATTATGGG	*	NH:i:1	HI:i:1	AS:i:29	nM:i:0
LEPR.11.0	0	chr1	65425374	255	5M140167N25M	*	0	0	TATGGGTGTACTTCTCTGAAGTAAGATGAT	*	NH:i:1	HI:i:1	AS:i:28	nM:i:0
LEPR.22.0	0	chr1	65525708	255	30M	*	0	0	GAGCGGCCCCATCGCAGAGCCCACGGCCAG	*	NH:i:1	HI:i:1	AS:i:29	nM:i:0
LEPR.19.0	0	chr1	65525775	255	30M	*	0	0	CGCCGCCATCTCTGCCTTCGGTCGAGTTGG	*	NH:i:1	HI:i:1	AS:i:29	nM:i:0
LEPR.17.0	0	chr1	65525795	255	24M39727N6M	*	0	0	GTCGAGTTGGACCCCCGGATCAAGGTGTAC	*	NH:i:1	HI:i:1	AS:i:28	nM:i:0
LEPR.11.0	0	chr1	65565544	255	30M	*	0	0	AGGTGTACTTCTCTGAAGTAAGATGATTTG	*	NH:i:1	HI:i:1	AS:i:29	nM:i:0
LEPR.10.0	0	chr1	65565576	255	30M	*	0	0	AAAAATTCTGTGTGGTTTTGTTACATTGGG	*	NH:i:1	HI:i:1	AS:i:29	nM:i:0
