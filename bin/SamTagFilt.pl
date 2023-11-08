#!/usr/bin/perl

$key="TP";
$limit=48;
($key)=@ARGV if @ARGV;
($key, $limit)=@ARGV if @ARGV>1;

while (<STDIN>) {
  if (/$key:i:(\d+)/) {
    next if $1>$limit;
  }
  print;
}
