#!/usr/bin/perl

foreach (@ARGV) {
  $G{$_}++;
}

$low='0'; #'!'
$high='z'; #'~'
$base=ord($low);
$limit=ord($high)-$base;
while (<STDIN>) {
  if (/>(\S+)/) { $gene=$1;
    $gene=$1 if /gene=(\S+)/;
    $flag=0;
    push(@lst, $seq) if $seq; undef $seq;
    if (exists $G{$gene}) {
      # print STDERR scalar(@lst)."\t$gene\n";
      s/^>/@/;
      push(@lst, $_);
      $flag=1;
      # check exons here
      die unless /exons:(\S+)/;
      $exons=$1;
      while ($exons=~/(\d+)-(\d+)/g) {
        ($s, $e)=($1, $2);
        foreach $n ($s .. $e) {
	  if ($CNT{$gene}{$n}>=$limit) {
            $LIMIT{$gene}++;
            warn "Gene $gene #exon depth reach limit[$limit]\n" if $LIMIT{$gene}==1;
          } else {
            $CNT{$gene}{$n}++;
          }
        }
      }
    }
  } elsif (/^[ACTGN]+$/i) { chomp();
    $seq.=uc($_) if $flag;
  } else {
    die $_;
  }
}
push(@lst, $seq) if $seq;

# print STDERR scalar(@lst);
#!(33)-~(126)=95
foreach (@lst) {
  # print STDERR $_; next;
  if (/^([ACTGN]+)$/i) { chomp(); $seq=$1;
    $nSeq=length($seq); $nQual=length($qual);
    die "[$nSeq:$nQual]" unless length($seq)==length($qual);
    print STDERR "$gene\t$tid\t".length($seq)."\n";
    print "$seq\n+\n$qual\n";
    next;
  }
  die $_ unless /^@(\S+)/; $tid=$1;
  die $_ unless /gene=(\S+)/; $gene=$1;
  die $_ unless /exons:(\S+)/; $exons=$1;
  die $_ unless /loc:\S+([+-])/; $strand=$1;
  print; 
  undef $qual;
  while ($exons=~/(\d+)-(\d+)/g) {
    ($s, $e)=($1, $2);
    foreach $n ($s .. $e) {
      $qual.=chr($base+$CNT{$gene}{$n});
    }
  }
  $qual=reverse $qual if $strand=~/-/;
}

__END__
>NR_046018 gene=DDX11L1 loc:chr1|11874-14409|+ exons:11874-12227,12613-12721,13221-14409 segs:1-354,355-463,464-1652
>NR_024540 gene=WASH7P loc:chr1|14362-29370|- exons:14362-14829,14970-15038,15796-15947,16607-16765,16858-17055,17233-17368,17606-17742,17915-18061,18268-18366,24738-24891,29321-29370 segs:1-50,51-204,205-303,304-450,451-587,588-723,724-921,922-1080,1081-1232,1233-1301,1302-1769
>NR_106918 gene=MIR6859-1 loc:chr1|17369-17436|- exons:17369-17436 segs:1-68
>NR_107062_3 gene=MIR6859-2 loc:chr1|17369-17436|- exons:17369-17436 segs:1-68
>NR_107063_1 gene=MIR6859-3 loc:chr1|17369-17436|- exons:17369-17436 segs:1-68
>NR_128720_1 gene=MIR6859-4 loc:chr1|17369-17436|- exons:17369-17436 segs:1-68
>NR_036051_3 gene=MIR1302-2 loc:chr1|30366-30503|+ exons:30366-30503 segs:1-138
>NR_036266_2 gene=MIR1302-9 loc:chr1|30366-30503|+ exons:30366-30503 segs:1-138
>NR_036267 gene=MIR1302-10 loc:chr1|30366-30503|+ exons:30366-30503 segs:1-138
>NR_036268_3 gene=MIR1302-11 loc:chr1|30366-30503|+ exons:30366-30503 segs:1-138
