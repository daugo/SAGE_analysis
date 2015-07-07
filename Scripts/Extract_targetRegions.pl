#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Bio::SeqIO;
use Bio::Seq;
use Bio::DB::Fasta;
use File::Basename;

my $supercontigs_fasta = ''; #Genomic sequences
my $transcripts_fasta = ''; #Predicted trascripts
my $annotation_file = ''; #Annotation file 
my $distance = ''; #Distance to extend
my $help; #help menu flag
my $license; #licencs information flag
my $version="0.01"; #Script version

#Get Options long argument definitions
GetOptions(
'genome_fasta|g=s' => \$supercontigs_fasta,
'transcripts_fasta|t=s' => \$transcripts_fasta,
'genome_feature_file|f=s' => \$annotation_file,
'distance|d=i' => \$distance,
'help|h' => \$help,
'license|l' => \$license, 
);

#Check help flag
if ($help) { 
	&usage();
	exit 0;
}
#Check license flag
if ($license) {
	&license();
	exit 0;
}
#Check that inputs files exist
unless (-s $supercontigs_fasta or -s $transcripts_fasta or -s  $annotation_file) {
	&usage();
	exit 0;
}
#Execute GetTargetSequences
&GetTargetSequences($supercontigs_fasta, $transcripts_fasta, $annotation_file, $distance);

sub GetTargetSequences {
	my $genome_fasta = shift;
	my $transcripts_fasta = shift;
	my $annotation_file = shift;
	my $distance = shift;
	
	my $hash_anno_ref;
	
	if ($annotation_file =~ /\.gff$/i) {
		$hash_anno_ref = &GetTranscriptPos_GFF($annotation_file);
	}
	elsif ($annotation_file =~ /\.gtf/i) {
		$hash_anno_ref = &GetTranscriptPos_GTF($annotation_file);
	}
	else {
		die "Anotation file name must have extension .gff or .gtf";
	}

	my ($fn, $dir, $suf) = fileparse($genome_fasta,qr/\.[^.]*/);
	my $outfile = $fn."SAGE_Targets"."_$distance".$suf;
	
	my $genome_fh = SeqIO_fasta_fh($genome_fasta,'<');
	
	my $transcript_db = Bio::DB::Fasta->new($transcripts_fasta);
	
	my $targets_seq = SeqIO_fasta_fh($outfile,'>');


	while (my $seq = $genome_fh->next_seq) {
		foreach my $transcript_id ( sort { $hash_anno_ref->{$seq->id}->{$a}{start} <=> $hash_anno_ref->{$seq->id}->{$b}{start} } keys %{$hash_anno_ref->{$seq->id}} ) {
			my $start = $hash_anno_ref->{$seq->id}->{$transcript_id}->{start};
			my $end = $hash_anno_ref->{$seq->id}->{$transcript_id}->{end};
			my $strand = $hash_anno_ref->{$seq->id}->{$transcript_id}->{strand};
			my $start_feature = $hash_anno_ref->{$seq->id}->{$transcript_id}->{feature_start};
			my $end_feature = $hash_anno_ref->{$seq->id}->{$transcript_id}->{feature_end};
			my $transcript_len = abs($end - $start);
			if ($end and $start) {
				my $subseq = '';
				my $extend_end = 0;
				#print "	$transcript_id\n";
				if ($strand eq '+') {
					if ($end + $distance  <= $seq->length) {
						$extend_end= $end + $distance;
					}
					else {
						$extend_end= $seq->length;
					}
					$subseq = $seq->subseq($end,$extend_end);
					$subseq = substr $subseq, 1-length($subseq);	
					
					my $transcript_obj = $transcript_db->get_Seq_by_id($transcript_id);
					my $transcript_seq = $transcript_obj->seq;
					$subseq = $transcript_seq . $subseq;
				}
				elsif ($strand eq '-') {
					if ($end - $distance > 0) {
						$extend_end = $end -$distance;
					}
					else {
						$extend_end = 1;
					}
					$subseq = $seq->subseq($extend_end,$end);
					$subseq = substr $subseq, 0, length($subseq)-1;
					my $transcript_obj = $transcript_db->get_Seq_by_id($transcript_id);
					my $transcript_seq = $transcript_obj->seq;
					my $transciptSeq = Bio::PrimarySeq->new(-seq => $transcript_seq,
										-id =>  $transcript_id,
										);
					$transciptSeq= $transciptSeq->revcom;
					$transcript_seq= $transciptSeq->seq;
					
					$subseq = $subseq.$transcript_seq ;
				}
				my $target_seq = Bio::PrimarySeq->new( -seq => $subseq,
												-id => "$transcript_id",
												-desc => "|transcript_len:$transcript_len|extend_len:$distance|3UTR($extend_end)|$start;$end;$extend_end|strand:$strand|start_feature:$start_feature|end_feature:$end_feature",
											);
				$targets_seq->write_seq($target_seq) if $strand eq '+';
				$targets_seq->write_seq($target_seq->revcom) if $strand eq '-';
			}
			else {
				warn "$transcript_id from superpercontig: ".$seq->id." doesn't have the start_codon and stop codon information";
			}
		}

	}
}


sub GetTranscriptPos_GTF {
	my $annotation_file = shift;
	open my $fh,'<', $annotation_file or die "Cannot open file $annotation_file";
	my %TrancriptPosHash;
	my %transcriptId_count;
	while (my $annon = <$fh>) {
		chomp($annon);
		my @fields = split('\t',$annon);
		my ($supercontig, $feature, $start, $end, $strand, $description) = ($fields[0], $fields[2], $fields[3], $fields[4], $fields[6], $fields[8]);
			my ($transcript_id) = $description =~ /; transcript_id "(\S+)";/;
		if ($strand eq '+') {
			if (!$transcriptId_count{$transcript_id} or $TrancriptPosHash{$supercontig}{$transcript_id}{'start'} > $start) {
				$TrancriptPosHash{$supercontig}{$transcript_id}{'start'} = $start;
				$TrancriptPosHash{$supercontig}{$transcript_id}{'feature_start'} = $feature;
			}
			if (!$transcriptId_count{$transcript_id} or $TrancriptPosHash{$supercontig}{$transcript_id}{'end'} < $end) {
				$TrancriptPosHash{$supercontig}{$transcript_id}{'end'} = $end;
				$TrancriptPosHash{$supercontig}{$transcript_id}{'feature_end'} = $feature;
			}
			$transcriptId_count{$transcript_id} += 1;
		}
		elsif ($strand eq '-') {
			if (!$transcriptId_count{$transcript_id} or $TrancriptPosHash{$supercontig}{$transcript_id}{'start'} < $end) {
				$TrancriptPosHash{$supercontig}{$transcript_id}{'start'} = $end;
				$TrancriptPosHash{$supercontig}{$transcript_id}{'feature_start'} = $feature;
			}
			if (!$transcriptId_count{$transcript_id} or $TrancriptPosHash{$supercontig}{$transcript_id}{'end'} > $start) {
				$TrancriptPosHash{$supercontig}{$transcript_id}{'end'} = $start;
				$TrancriptPosHash{$supercontig}{$transcript_id}{'feature_end'} = $feature;
			}
			$transcriptId_count{$transcript_id} += 1;
		}
		$TrancriptPosHash{$supercontig}{$transcript_id}{'strand'} = $strand;
	}
	return \%TrancriptPosHash;
}


sub GetTranscriptPos_GFF {
	my $annotation_file = shift;
	open my $fh,'<', $annotation_file or die "Cannot open file $annotation_file";
	my %TrancriptPosHash;
	#my %transcriptId_count;
	while (my $annon = <$fh>) {
		chomp($annon);
		my @fields = split('\t',$annon);
		my ($supercontig, $feature, $start, $end, $strand, $description) = ($fields[0], $fields[2], $fields[3], $fields[4], $fields[6], $fields[8]);
		my ($transcript_id) = $description =~ /ID=([0-9A-Z]+);/i;
		
		if ($feature eq 'mRNA') {
			if ($strand eq '+') {
				$TrancriptPosHash{$supercontig}{$transcript_id}{'start'} = $start;
				$TrancriptPosHash{$supercontig}{$transcript_id}{'feature_start'} = $feature;
				$TrancriptPosHash{$supercontig}{$transcript_id}{'end'} = $end;
				$TrancriptPosHash{$supercontig}{$transcript_id}{'feature_end'} = $feature;
			}
			elsif ($strand eq '-') {
				$TrancriptPosHash{$supercontig}{$transcript_id}{'start'} = $end;
				$TrancriptPosHash{$supercontig}{$transcript_id}{'feature_start'} = $feature;
				$TrancriptPosHash{$supercontig}{$transcript_id}{'end'} = $start;
				$TrancriptPosHash{$supercontig}{$transcript_id}{'feature_end'} = $feature;
			}
			$TrancriptPosHash{$supercontig}{$transcript_id}{'strand'} = $strand;
		}
	}
	return \%TrancriptPosHash;
}


sub SeqIO_fasta_fh {
		my $fasta_file = shift;
		my $io = shift;
		die "Not valid input/ouput sign for Bio::SeqIO->new" if $io !~ /^[><]$/ ;
		my $fh = Bio::SeqIO->new(
						-format => 'fasta',
						-file => $io.$fasta_file,
						);
		return $fh;
}


sub usage{
    print STDERR "$0 version $version, Copyright (C) 2013 David Urbina-Gomez\n";
    print STDERR "$0 comes with ABSOLUTELY NO WARRANTY; for details type `$0 -l'.\n";
    print STDERR "This is free software, and you are welcome to redistribute it under certain conditions;\n";
    print STDERR "type `$0 -l' for details.\n";
    print STDERR <<EOF;
NAME
    Extract_targetRegions.pl  Extend predicted transcrits to a certain length at the 3' end (needed a gff or gtf file and a fasta of the genomic regions)

USAGE
    Extract_targetRegion.pl -g [genome.fasta] -t [transcripts.fasta]-f [genome features file(gff|gtf)] -d [maximum distance to extend ]
 
OPTIONS
    --genome_fasta,         -g    Genome fasta                              REQUIRED
    --transcripts_fasta,    -t    Transcripts fasta                         REQUIRED
    --genome_feature_file,  -f    Genome feature file (gtf|gff)             REQUIRED
    --distance,             -d    maximum distance to extend the 3' end     REQUIRED
    --help,                 -h    This help
    --license,              -l    License.

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
