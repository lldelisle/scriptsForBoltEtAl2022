#!/usr/local/bin/perl -w
# NOTE: YOU MUST CHANGE THE LINE ABOVE TO POINT TO
# THE FULL PATH OF THE PERL EXECUTABLE ON YOUR SYSTEM.

# Please see copyright notice and system requirements
# in this document. 

# This program used to produce rectangular dot plots for:
# 
# Hughes, J.F., Skaletsky, H., Pyntikova, T., Minx, P.J.,
# Graves, T., Rozen, S., Wilson, R.K., Page, D.C.
# Conservation of Y-linked genes during human evolution
# revealed by comparative sequencing in chimpanzee.
# Nature 437, 100-103 (2005)

# Copyright (c) 2008 Helen Skaletsky and Whitehead Institute
# 
# Redistribution and use in source and binary forms, with or
# without modification, are permitted provided that the
# following conditions are met:
#
# 1.  Redistributions must reproduce the above copyright
# notice, this list of conditions and the following
# disclaimer in the documentation and/or other materials
# provided with the distribution.  Redistributions of source
# code must also reproduce this information in the source code
# itself.
#
# 2.  If the program is modified, redistributions must
# include a notice (in the same places as above) indicating
# that the redistributed program is not identical to the
# version distributed by Whitehead Institute.
#
# 3.  All advertising materials mentioning features or use
# of this software must display the following
# acknowledgment:
#
#         This product includes software developed by 
#         Helen Skaletsky and the Whitehead Institute
#         for Biomedical Research.
#
# 4.  The name of the Whitehead Institute may not be used to
# endorse or promote products derived from this software
# without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE WHITEHEAD INSTITUTE AND
# HELEN SKALETSKY ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE WHITEHEAD
# INSTITUTE OR HELEN SKALETSKY BE LIABLE FOR ANY DIRECT,
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
# IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# System requirements: requires perl and perl module GD.pm,
# available from CPAN (author LDS --- Lincoln D. Stein).
 
# For usage information, run this program with the flag
# -h.

require v5.6.0;

use GD;
use strict;
use vars qw/$VERSION/;

BEGIN {

    $VERSION = '1.0';

    if (!defined $GD::VERSION || $GD::VERSION ne 1.32) {
	print STDERR "\n$0:\n";
	print STDERR "WARNING -- GD VERSION 1.32 REQUIRED.\n";
	print STDERR "You have GD version $GD::VERSION.\n"
	    if defined $GD::VERSION;
	print STDERR "If $0\n";
	print STDERR "does not work correctly you will\n";
	print STDERR "have to install GD version 1.32.\n\n";
    }
}

use Getopt::Long;

sub usage();
sub main();
sub print_dots ($$$$$$$$$);
sub print_png ($$);

main();

sub usage() {
    print qq/
USAGE: $0 -w <word size> -s <step length>
          -1 <1st seq> -2 <2nd seq> -o <out file> 
          -d <dot file> -t <title> [ -h ]

Create a PNG ("portable network graphics") file
that displays a rectangular dot plot of the first 
input sequence against the second input sequence.

<word size> is the word size for a match.  A dot
            is printed if there is a perfect match
            of length <word size>.

<step length> is the number of bases to move the
            word for each dot.

<in file 1>   is a fasta format file from which the
            1st sequence is taken.

<in file 2>   is a fasta format file from which the
            2nd sequence is taken.

<out file>  is the PNG file created.

<dot file>  contains 0 based positions of perfect
            matches of length <word size>.
            (E.g., the line "251 1077" means that
            substrings of length <word size>
            starting at 251 in the 1st sequence and 
            at 1077 in the 2nd sequence are
            identical.)

<title>     is a title to place in the output.

-h causes this message to printed.

(Version $VERSION)

/;
}

sub main() {

    my ($seqfile1, $seqfile2, $word, $step, $outfile, $title, $dotfile, $help);

    if (!GetOptions('1st=s'  => \$seqfile1,
                    '2nd=s'  => \$seqfile2,
		    'wordlen=i' => \$word,
		    'step=i'    => \$step,
                    'outfile=s' => \$outfile,
                    'title=s'   => \$title,
		    'dotfile=s' => \$dotfile,
		    'help'      => \$help)) {
	usage;
	exit -1;
    }
    if ($help) {
	usage;
	exit;
    }

    if (!defined $word
	|| !defined $step
	|| !defined $seqfile1
        || !defined $seqfile2
	|| !defined $outfile
	|| !defined $dotfile) {
	usage; exit -1;
    }
	    
    $title = '' unless defined $title;

    if ($word <= 0) {
	print STDERR "$0 -w <word size> must be >= 0\n";
	exit -1;
    }

    if ($step <= 0) {
	print STDERR "$0 -s <step length> must be >= 0\n";
	exit -1;
    }


    open(IN, $seqfile1) || die "Cannot open $seqfile1: $!";
    <IN>;
    my $seq1 = '';
    while(<IN>){
	chop;  $seq1 .= uc($_);
    }
    close IN;
    my $m1 = length($seq1);

    open(IN, $seqfile2) || die "Cannot open $seqfile2: $!";
    <IN>;
    my $seq2 = '';
    while(<IN>){
        chop;  $seq2 .= uc($_);
    }
    close IN;
    my $m2 = length($seq2);

    my $k = 0; 
    my $n = $m1 - $word;

    open(OUT, ">$dotfile")
	|| die "Cannot write $dotfile: $!\n";

    while($k < $n) {
	my $s = substr($seq1, $k, $word);
        if($s =~ /N/){
            $k += $step;
            next;
        }
	while($seq2 =~ m/$s/g){
	    my $t = pos $seq2; $t -= $word;
	    print OUT "$k\t$t\n";
	}
	my $s1 = reverse($s);
	$s1 =~ tr/ACGTURYMWSKDHVBN/TGCAAYRKWSMHDBVN/;
	while($seq2 =~ m/$s1/g){
	    my $t = pos $seq2; $t -= $word;
	    print OUT "$k\t$t\n";
	}
	$k += $step;
    }
    close OUT;

    # Create and print the output.
    my $width = 700; my $height = 730; my $x0 = 30;
    my $img = new GD::Image($width, $height);
    my $white = $img->colorAllocate(255,255,255);
    my $black = $img->colorAllocate(0,0,0);
    my $gray = $img->colorAllocate(187,187,187);
    $img->interlaced('true');

    $x0 = 60; 
    $img->string(gdLargeFont, $x0, 30, 
		 "$title -w=$word -s=$step",
		 $black);
    print_dots($img, $width, $height, $x0, $black, $gray, $m1, $m2, $dotfile);
    print_png($img, $outfile)
}

sub print_dots ($$$$$$$$$) {
    my ($img, $width, $height, $x0, $black, $gray, $m1, $m2, $dotfile) = @_;
    open(IN, $dotfile) || die "Cannot open $dotfile: $!\n";
    my ($x1, $y1, $del, $m);
    if($m1 > $m2){$m = $m1;}
    else {$m = $m2;}
    $del = ($width - 80)/$m;
    my $l = $del * $m1;
    my $h = $del * $m2;
    my $count = 0;
    while(<IN>){
	chop; 
	my @x = split '\t';
	# For progress report on very long run
        # if(int($count/10000)*10000 == $count){
        #    print STDERR "$x[0]   $x[1]\n";
	# }
	$x1 = $x0 + $del*$x[0];
	$y1 = $x0 + $del*$x[1];
	$img->setPixel($x1, $y1, $black);
	$count++;
    }
    close IN;
    $img->line($x0 + $l, $x0, $x0 + $l, $x0 + $h, $gray);
    $img->line($x0, $x0 + $h, $x0 + $l, $x0 + $h, $gray);
    $img->line($x0, $x0, $x0, $x0 + $h, $gray);
    $img->line($x0, $x0, $x0 + $l, $x0, $gray);
    $img->string(gdLargeFont, $x0, $x0 + $h + 5,
                 "$m1 bp", $black);
    my $b;
    if($x0 < 18){
        $b = 0;
    }
    else{
        $b = $x0 - 18;
    }
    $img->stringUp(gdLargeFont, $b, $x0 + $h, 
                 "$m2 bp", $black);
}

sub print_png ($$) {
    my ($img, $outfile) = @_;
    if ($outfile) {
	open(OUT, ">$outfile")
	    || die "Cannot write $outfile: $!\n";
	print OUT $img->png;
	close OUT;
    } else {
	print $img->png;
    }
}


