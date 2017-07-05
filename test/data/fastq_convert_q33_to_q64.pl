#!/usr/bin/perl


#def qual33(qual64): return chr(ord(qual64)-31)
#def qual64(qual33): return chr(ord(qual33)+31)


while (<STDIN>) {
    chomp;
    if ($_ =~ /^\@/) {
	print $_."\n";
        $_ = <>;
	chomp;
	print $_."\n";
	$_ = <>;
	chomp;
	print $_."\n";
	$qual33_line = <>;
	chomp $qual33_line;
	$qual64_line = '';
	foreach $q33 (split (//, $qual33_line)) {
	    $q64 = chr(ord($q33)+31);
	    $qual64_line .= $q64;
	}
	print $qual64_line."\n";
    }
}
