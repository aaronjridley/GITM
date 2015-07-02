#!/usr/bin/perl
#  Copyright (C) 2002 Regents of the University of Michigan, portions used with permission 
#  For more information, see http://csem.engin.umich.edu/tools/swmf

# (c) Gabor Toth 2007

sub XmlRead{

    # The XmlRead function implements an extremely light weight XML reader.
    # It is a single Perl subroutine that does not use any package or C code.
    #
    # Input is a string with XML text and the function returns a Perl data
    # structure with the parsed XML text. The sturcture of the returned value
    # follows the XML::Parser::EasyTree package.
    #
    # The returned value is a reference to an array of hashes. 
    # Each hash has a key 'type' with value either 't' (text) or 'e' 
    # (XML element). For type 't' the hash has another key 'content'
    # with a value containing the plain text.
    # For type 'e' the key 'name' contains the name of the tag.
    # If the tag has attributes the key 'attrib' has a value that is a 
    # reference to a hash of attributes. The attribute hash has keys that
    # are the attribute names and values that are the attribute values.
    # If the opening tag is paired with a closing tag, the hash has a 
    # key 'content' with a value that is a reference to an array of hashes that
    # describes the XML text between the opening and closing tags
    # as defined above.
    #
    # Although this is a recursive definition, the returned tree 
    # always has a finite depth since the input text is of finite length.

    use strict;
    use vars '@Xml'; # because of 'local' that is needed for recursiveness

    local($_) = shift;  # XML text
    local(@Xml);        # Array of XML elements
    my $KeepSpace = shift;

    # Debugging flag
    my $Debug = 0;

    # XML special characters
    my %unxml = ( 'gt' => '>',
		  'lt' => '<',
		  'amp'=> '&',
		  'quot' => '"',
		  'apos' => "'",
		  );

    my $ERROR = "ERROR in XmlRead.pl";

    # Remove all comments
    s/<!--.*?-->//gs;

    # Find XML elements <...>
    while( /^([^<]*)</ ){
	my $text = $1; $_ = $';

	$text =~ s/^\s*$// unless $KeepSpace;

	# Save text before next <
	if($text){
	    $text =~ s/\&(\w+);/$unxml{$1}/g;
	    push(@Xml, { 'type' => 't', 'content' => $text } );
	    warn "text = $text\n" if $Debug;
	}

	# Now read the XML tag
	my $ElementRef;
	if( /^(\w+)\s*/ ){
	    my $tag = $1; $_ = $';
	    $ElementRef->{type} = "e";
	    $ElementRef->{name} = $tag;

	    # Read attributes
	    while( /^(\w+)\s*=\s*(\'[^\']*\'|\"[^\"]*\")\s*/ ){
		my $attribute = $1; $_ = $';
		my $value = $2; $value =~ s/^[\'\"]//; $value =~ s/[\'\"]$//;
		$value =~ s/\&(\w+);/$unxml{$1}/g;

		$ElementRef->{attrib}{$attribute} = $value;

		warn "attribute=$attribute, value=$value\n" if $Debug;
	    }

	    # Check if the last character is a /
	    my $endtag = (s|^/||);

	    # Check if the closing > is there
	    if(not s/^>//){
		warn "$ERROR: incorrect attribute for <$tag ... ".
		    substr($_,0,100)."\n";
		return;
	    }

	    warn "tag=$tag, endtag=$endtag\n" if $Debug;

	    if(not $endtag){
		# Check for nested <tag...> of the same type and </tag>
		my $content;
		my $level = 1; # there was an opening tag so far
		while( /(<$tag\b[^>]*>|<\/$tag>)/ ){
		    $content .= $`; $_ = $';
		    my $tag2 = $1;
		    if($tag2 eq "</$tag>"){
			$level--;
			last unless $level; # exit if we got back to level 0
		    }elsif($tag2 !~ /\/>$/){
			$level++;
		    }
		    $content .= $tag2; # still inside, keep searching for end
		}
		if($level){
		    warn "$ERROR: missing </$tag> after <$tag ...>".
			substr($content.$_,0,100)."\n";
		    return;
		}

		$content =~ s/^\s*$// unless $KeepSpace;

		$ElementRef->{content} = &XmlRead($content, $KeepSpace) 
		    if length $content;
	    }
	}else{
	    warn "$ERROR: incorrect tag at <".substr($_,0,100)."\n";
	    return;
	}
	push(@Xml, $ElementRef);
    }

    s/^\s*$// unless $KeepSpace;

    # Check for closing text
    if( not /^\n?$/ ){
	s/\&(\w+);/$unxml{$1}/g;
	push(@Xml, { 'type' => 't', 'content' => $_ } );
	warn "final text = $_\n" if $Debug;
    }

    # Return a reference to the @Xml array
    return [ @Xml ]; 
}

1;
