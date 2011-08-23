#!/homes/gws/aritter/local/bin/perl -w

use strict;

use Getopt::Long;
use List::Util qw[max min shuffle];

#Need at least 50 tuples per relation
use constant MIN_TUPLES => 20;

my %arg1docs;
my %arg2docs;

my %arg1vocab;
my %arg2vocab;

my $a1next = 1;	#NOTE: it is very important that these start at 1 (and not 0!)
my $a2next = 1;

my ($tupleFile, $a1vocab, $a2vocab, $maxVocSize, $maxTuplesPerRel);
GetOptions(
    'tupleFile=s'  => \$tupleFile,
    'a1vocab=s'    => \$a1vocab,
    'a2vocab=s'    => \$a2vocab,
    'maxTuplesPerRel=i' => \$maxTuplesPerRel,
    'maxVocSize=i' => \$maxVocSize
    );

if(!$tupleFile || ($a1vocab && !$a2vocab) || (!$a1vocab && $a2vocab)) {
    die "usage: $0 --tupleFile=X [--a1vocab=X --a2vocab=X] [--maxTuplesPerRel=X] [--maxVocSize=X]\n";
}

open TUPLES, $tupleFile;

#############################################################################
# Predefined vocabulary 
#############################################################################
if($a1vocab && $a2vocab) {
    open A1_VOCAB, $a1vocab;
    open A2_VOCAB, $a2vocab;

    while(<A1_VOCAB>) {
        chomp;
        my ($string, $id) = split /\t/;
        $arg1vocab{$string} = $id;
        if ($id >= $a1next) {
            $a1next = $id + 1;
        }
    }

    while(<A2_VOCAB>) {
        chomp;
        my ($string, $id) = split /\t/;
        $arg2vocab{$string} = $id;
        if ($id >= $a2next) {
            $a2next = $id + 1;
        }
    }

    #This should cause an error if there is something outside the vocab
    #$a1next = -100;
    #$a2next = -100;

    close A1_VOCAB;
    close A2_VOCAB;
}

sub arg1v {
    my $string = shift;
    if(!$arg1vocab{$string}) {
        if($maxVocSize && $a1next == $maxVocSize) {
            #Once the vocab is full, everything else is OOV
            $string = "-OOV-";
            $arg1vocab{$string} = $a1next;
        } else {
            $arg1vocab{$string} = $a1next;
            $a1next++;
        }
    }
    return $arg1vocab{$string};
}

sub arg2v {
    my $string = shift;
    if(!$arg2vocab{$string}) {
        if($maxVocSize && $a2next == $maxVocSize) {
            #Once the vocab is full, everything else is OOV
            $string = "-OOV-";
            $arg2vocab{$string} = $a2next;
        } else {
            $arg2vocab{$string} = $a2next;
            $a2next++;
        }
    }
    return $arg2vocab{$string};
}

while(<TUPLES>) {
    chomp;
    my ($arg1, $pred, $arg2) = split /\t/;
    push @{$arg1docs{$pred}}, arg1v($arg1);
    push @{$arg2docs{$pred}}, arg2v($arg2);
}

open PREDS, ">preds";
open ARG1, ">arg1";
open ARG2, ">arg2";

my @preds = keys %arg1docs;

my @a1;
my @a2;
for my $pred (@preds) {
    next if(@{$arg1docs{$pred}} < MIN_TUPLES);
    print PREDS "$pred\n";
    @a1 = shuffle @{$arg1docs{$pred}};
    @a2 = shuffle @{$arg2docs{$pred}};
    if($maxTuplesPerRel) {
        print ARG1 join(" ", @a1[0..min($maxTuplesPerRel, $#a1)]) . "\n";
        print ARG2 join(" ", @a2[0..min($maxTuplesPerRel, $#a2)]) . "\n";
    } else {
        print ARG1 join(" ", @a1) . "\n";
        print ARG2 join(" ", @a2) . "\n";
    }
}

open A1VOCAB, ">a1vocab";
for my $word (keys %arg1vocab) {
    print A1VOCAB "$word\t$arg1vocab{$word}\n";
}

open A2VOCAB, ">a2vocab";
for my $word (keys %arg2vocab) {
    print A2VOCAB "$word\t$arg2vocab{$word}\n";
}
