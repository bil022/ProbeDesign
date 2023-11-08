#!/usr/bin/perl

die "Usage: $0 <usr_lst> <ref_genes>" unless @ARGV;

($usr_lst, $genes)=@ARGV;
open GENES, "<$genes" or die "$genes?";
while (<GENES>) {
  ($id)=split();
  $GN{uc($id)}=$id;
}

open USR, "<$usr_lst" or die "gene_list.txt?[$usr_lst]";
while (<USR>) {
  ($id)=split(); $total++;
  if (exists $GN{uc($id)}) {
    print $GN{uc($id)}."\n";
    next;
  }
  warn "$id not found";
  $warn++;
}

if ($warn) {
  print STDERR "$warn/$total genes not found\n";
}

__END__
GFAP
S100b
ALDH1
LCN2
Glast
Gli1
Hes5
HopX
Vim
Steap4
S1pr3
Serpin3n
Serping1
PTBP1
Aqp4
Serpinf1
Gja1
Mfge8
Slco1c1
Aldoc1
GLT1
AGT
Mertk
Unc13c
MCM2
Ki67
Tbr2
ASCL1
Oct3/4.
Kif4
c-Myc
SOX2
Nestin
Fam107a
Tbr1
NFIA
SOX11
SOX9
GFAP
Glast
BLBP
Pax6
Vimentin
ASCL1
Nestin
Hes1
Hes5
Frzb
Slc1a3
Mfge8
mertk
DCX
CRMP4
NeuroD1
DLX2
LHX6
Pax6
TBR1
Prox
NeuN
PTBP2
Glut1
GABA
miR9
miR124
Vip
Sst
Pvalb
Thy1
Gad1
Spink8
Gm11549
Pnoc
Pcp4
Neurogn2
Neurogn1
DLX1
Arx
Calb1
FoxA1
FoxA2
Gpc4
Vps13c
Wfs1
Map3k15
Pvrl3
Dsp
Slit2
Nmb
Grp
Hapln2
Olig1
Olig2
MBP
Sox10
pdgfra
IBA1
Mgl1
Mgl2
Pvm1
Pvm2
Pmac
Aif1
MRF1
IBA2
cldn5
Ly6c1
Acta2
Aldoc
Smarca4
Sin3a
Npas3
Neurod4
tdTomato
TDP-43
HTT
Malat1
MyoD
Myf5
