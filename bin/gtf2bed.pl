#!/usr/bin/perl

if (@ARGV) {
  while (<>) {
    next unless /\texon\t/;
    ($chr, $src, $exon, $s, $e, $a, $strand, $b, $g, $gid)=split();
    die unless $g eq "gene_id";
    $gid=~s/[\";]//g;
    die $_ unless /gene_name "(\S+)";/; $gid=$1;
    print "$chr:+:$gid\t$s\t$e\n" if $strand=~/\+/;
    print "$chr:-:$gid\t$e\t$s\n" if $strand=~/-/;
  }
} else {
  die "hg38?";
  while (<>) {
    die unless /(\S+):(\S):(\S+)\t(\d+)\t(\d+)/;
    ($chr, $strand, $gid, $s, $e)=($1,$2,$3,$4,$5,$6);
    $N++; $N=1 if $gid ne $prev;
    print "\n" if $gid ne $prev;
    #print "$chr\tunknown\texon\t$s\t$e\t.\t$strand\t.\tgene_id \"$gid\";\n";
    printf "samtools faidx hg38.fa" if $N==1;
    printf " $chr:$s-$e";
    $prev=$gid;
  }
  print "\n" if $gid;
}

__END__
awk '{print $1$7$10"\t"$4"\t"$5}' exons.gtf | sed 's/[\";]//g' | sort -k1,1 -k2,2n | bedtools merge -i - ^C
[bli@silencer hg38]$ head exons.gtf 
chr19	unknown	exon	58346806	58347029	.	-	.	gene_id "A1BG"; gene_name "A1BG"; p_id "P8460"; transcript_id "NM_130786"; tss_id "TSS7852";
