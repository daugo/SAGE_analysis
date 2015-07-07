#!/usr/bin/perl
#SAGE_inSilico.pl

use strict;
use warnings;
use Getopt::Long;

use Bio::Restriction::Enzyme;
use Bio::Restriction::Analysis;
use Bio::SeqIO;
###############################
#Get option from command line
###############################


my $version='0.01';
my $license='';
my $help='';
 

my $seqdb;
my $outfasta='';
my $enzyme_name;
my $recognition;
my $cut_site='';
my $sage_fragment_size;


GetOptions(
    'license|l'      => \$license,
    'help|h|?'       => \$help,
    'in_file|i=s'      => \$seqdb,
    'outfile|o=s' => \$outfasta,
    'enzyme_name|e=s'   => \$enzyme_name,
    'recognition_seq|r=s' => \$recognition,
    'cut_site|c=s' => \$cut_site,
    'fragment_size|f=i' => \$sage_fragment_size,

);

unless (-s $seqdb ) {
	&usage();
	die "$seqdb not found";
}

if (-s $outfasta) {
	warn "$outfasta exists: Overwriting";
}

&usage() if ($help);
&licence() if ($license);

restrictionAnalysis($seqdb, $outfasta, $enzyme_name, $recognition, $cut_site, $sage_fragment_size);

sub restrictionAnalysis {
	
	my $seqdb = shift;
	my $outfasta = shift;
	my $enzyme_name = shift;
	my $recognition =  shift;
	my $cut_site = shift;
	my $sage_fragment_size = shift;
	
	
	#Define Restriction Enzyme Object
	my $re = Bio::Restriction::Enzyme->new(-enzyme=>$enzyme_name,-seq=>$recognition,-cut=>$cut_site);

	#SeqIO filehandle (Read) 
	my $fasta_in_fh = Bio::SeqIO->new(
			-format => 'fasta',
			-file => "<$seqdb",
			);
	#SeqIO filehandle (Write)
	my $fasta_out_fh = Bio::SeqIO->new(
			-format => 'fasta',
			-file => ">$outfasta",
			);
			  
	while (my $seq = $fasta_in_fh->next_seq) {
		my $ra = Bio::Restriction::Analysis->new(-seq=>$seq, -enzymes=>$re);
		my @fragments = $ra->fragments($enzyme_name);
		
		for (0..$#fragments) {
			#revisar print secuencias en Extract_target
			my $fragment_seq = Bio::PrimarySeq->new( -seq => substr($fragments[$_],0, $sage_fragment_size),
							-id => $seq->id.";$_",
							-desc => $seq->desc,
							);
			$fasta_out_fh->write_seq($fragment_seq);
		}
	}

}


sub usage{
    print STDERR "$0 version $version, Copyright (C) 2013 David Urbina-Gomez\n";
    print STDERR "$0 comes with ABSOLUTELY NO WARRANTY; for details type `$0 -l'.\n";
    print STDERR "This is free software, and you are welcome to redistribute it under certain conditions;\n";
    print STDERR "type `$0 -l' for details.\n";
    print STDERR <<EOF;
NAME
    $0  Perform and in silico SAGE with certain enzyme Type II and a size of fragment obtained

USAGE
    $0 -i [in_file] -o [outfile] -e [enzyme_name] -c [cut_site] -r [recognition_seq] -f [fragment_size] 
 
OPTIONS
    --in_file,         -i    Transcripts fasta sequences                    REQUIRED
    --outfile,         -o    Sage possible fragments for each transcript    REQUIRED
    --enzyme_name,     -e    Enzyme name				    REQUIRED
    --recognition_seq  -r    Recognition sequence			    REQUIRED
    --cut_site         -c    Enzyme cut site				    OPTIONAL*
    --fragment_size    -f    SAGE fragment size				    REQUIRED
    --help,            -h    This help										
    --license,         -l    License.

EOF
    exit;
}
	
sub license{
    print STDERR <<EOF;

Copyright (C) 2013 David Urbina-Gomez
e-mail: daugo182\@gmail.com

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
EOF
exit;
}

exit 0;

