#!/usr/bin/perl -w

sub max($$) {if($_[0] > $_[1]) {return $_[0];} return $_[1]; }
sub min($$) {if($_[0] < $_[1]) {return $_[0];} return $_[1]; }

sub log10 {my $n = shift; return log($n)/log(10);}

sub readpqrconf(){

  #check for EOF

  if (eof(DAT)){
    return 0;
  }

  ++$nconf;
  @conf = (); # initialize config

  # read headers ie the REMARK lines
  for($i = 0;$i < $NREMARK; ++$i){
    $conf[$i] = <DAT>;
    if(!($conf[$i] =~ m/REMARK/)){
      die "Something wrong with header in conf # $nconf line # $i\n";
    }
  }

  if($conf[1] =~ m/total_atoms=([0-9]+)/){ # parse out natoms
    $natoms = $1;
  } else {
    die "Something when wrong while reading in header in conf # $nconf \n"
  }
  # read in atoms lines
  for($i = 0;$i < $natoms ; ++$i){
    $conf[$i+$NREMARK] = <DAT>;
    if(!($conf[$i+$NREMARK] =~ m/ATOM/)){
      die "Something wrong with ATOMs in conf # $nconf (atom # $i)\n";
    }
  }

  # read headers ie the REMARK lines
  for($i = 0;$i < $MREMARK; ++$i){
    $conf[$i+$NREMARK+$natoms] = <DAT>;
    if(!($conf[$i+$NREMARK+$natoms] =~ m/REMARK/)){
      die "Something wrong with header in conf # $nconf line # $i\n";
    }
  }

  # read in endmdl
  $conf[$natoms+$NREMARK+$MREMARK] = <DAT>;
  if(!($conf[$natoms+$NREMARK+$MREMARK] =~ m/^ENDMDL/)){
      die "Something wrong with ENDMDLs in conf # $nconf)\n";
    }
  return 1;
}

sub writepqrconf(){
  if ($nconf == 1){
    open (OUTPUT, ">$outfile") || die ("Could not open the output file $outfile!");
    print "\ncreating new file $outfile\n";
  }
  $n = 0;
  foreach $line (@conf){
    if(!$line){
      print "Blank line at $n?\n";
    }
    print OUTPUT $line;
    ++$n;
  }
}

# ----------------------- Begin Main ----------------------------------

# defaults

$NREMARK=4; # number of remark lines at beginning (could change to count them)
$MREMARK=3; # number of remark lines at end of file
$outfile = "output.pqr";
$nconf = 0;
$nconf1 = 0;
$data_file = $ARGV[0];
$!; # flush stdout

if($ARGV[1]){ #set output if necessary
  $outfile = $ARGV[1];
  print "set output to $outfile\n";
}

print "parsing $data_file\n";

open (DAT, $data_file) || die ("Could not open the file!");

$nmax = -1;
$nmin = -1;

print "Initial loop of $data_file\n";
# loop over file to get max atoms
while ($line = <DAT>) {
  if($line =~ m/REMARK/){
#    print $line;
    if($line =~ m/total_atoms=([0-9]+)/){
      if($nmax == -1){ $nmax = $nmin = $1;}
      $nmax = max($1,$nmax);
      $nmin = min($1,$nmin);
      ++$nconf1;
    }
  }
}

print "maximum atom = $nmax, minimum atoms = $nmin, nconf = $nconf1\n";

#rewind file
seek(DAT, 0, 0) or die "Can't seek to beginning of file: $!";

while (readpqrconf()){

  if($conf[1] =~ m/total_atoms=([0-9]+)/){ # parse out natoms
    $natoms = $1;
    chomp $conf[1];
    $conf[1] = $conf[1].", now=$nmax\n";
  } else {
    die "Something went wrong while reading in header in conf # $nconf \n"
  }
  print "\rnconf = $nconf/$nconf1, max = $nmax, natoms = $natoms         ";
  #move endline and other REMARKs to end
  for($i = 0;$i < $MREMARK+1; ++$i){
    $conf[$nmax+$NREMARK+$i] = $conf[$natoms+$NREMARK+$i];
  }

  # create new line to append - 
  # same x,y,z but molecule numer while preserving space.
  $line = $conf[$natoms-1+$NREMARK]; # initial line
  @field=split(/\s+/,$line); #split on white spaces

  $nmol = $field[5]+1; # get molecule number
  $line =~ m/(\s$field[4]\s+)$field[5]/; # pull out field with spaces
  $anew = $1;
  $anew =~ s/M/A/; # replace M with A
  # check for added characture
  if(length("$nmol") != length("$field[5]")){
    $anew =~ s/ //; # delete space if necessary
  }
  # create new line to copy
  $line =~ s/(\s$field[4]\s+)$field[5]/$anew$nmol/; #sub in the new number

  for($i=$natoms+$NREMARK;$i<$nmax+$NREMARK;++$i){
    $conf[$i] = $line; # copy atom line to fill out array
  }
  writepqrconf();
}

#clean up files
close DAT;
close OUTPUT;

print "\nWrote $nconf configurations\n";

exit(0);
