
###perl parser to get gene of unknown mgi

$mgi = $ARGV[0];
#$mgi='MGI:5613649';

my $url=`curl -s "http://www.informatics.jax.org/allele/$mgi"`;

if( $url=~ /Gene:<\/font(\s|.)+class='MP'>(.+)<\/a>&nbsp;&nbsp;/ ) {
	print STDOUT $2;
} else {
	print STDOUT "no_match";
}

