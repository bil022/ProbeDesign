#!/usr/bin/perl

while (<DATA>) {
  next if /^#/;
  ($id, $limit)=split();
  $LIMIT{uc($id)}=$limit;
}

while (<>) {
  if (/^@/) { print; next; }
  die unless /^(\S+)\.\d+\.\d+\.\d+\.\d+\.\d+\s/; $id=uc($1);
  die unless /PW:i:(\d+)/; $pw=$1;
  $limit=48; $limit=$LIMIT{$id} if exists $LIMIT{$id};
  print if $pw<=$limit;
}

__END__
#Malat1	60
#Htt	60
#Ptbp1	60
