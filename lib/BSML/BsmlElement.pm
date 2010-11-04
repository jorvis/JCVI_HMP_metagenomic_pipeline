package BSML::BsmlElement;

use warnings;
use strict;
use Data::Dumper;

BEGIN {
use BSML::Logger;
}


my $logger = BSML::Logger::get_logger("Logger::BSML");


=head1 NAME

BsmlElement.pm - Base Class for Bsml Elements

=head1 VERSION

This document refers to version 1.0 of the BSML Object Layer

=head1 SYNOPSIS

  $bsmlelem->setattr( 'id', 'A11243' );
  $hashref = $bsmlelem->returnattrHashR();

  foreach my $key ( keys( %{$hashref} ))
  {
    print "$hashref->{$key}\n";
  }

=head1 DESCRIPTION

=head2 Overview

  This file provides a base class for handling element attributes and BSML elements.

=head2 Constructor and initialization 

  A BsmlElement is not useful by itself. It should be only used as a base class.

=head2 Class and object methods

=over 4

=cut

=item $elem->addattr( $key, $value )

B<Description:> add an attribute to a Bsml element

B<Parameters:> ($key, $value) - a key value pair 

B<Returns:> None

=cut

sub addattr
  {

      $logger->debug("") if $logger->is_debug;


    my $self = shift;
    my ( $key, $value ) = @_;
    if( defined( $value ) ){

        ## clean up id attributes when adding
        if ($key eq 'id') {
           my $clean_id = getCleanID($value);
           if ($clean_id ne $value) {
                if (   $self->isa('BSML::BsmlFeature')
                    || $self->isa('BSML::BsmlSequence')) {
                    push(@{$self->{ 'BsmlAttr' }->{ 'primary_label' }}, $value);
                }
               $value = $clean_id;
           }
        ## clean up attributes which reference ids
        } elsif (
                $key eq 'refseq' 
            ||  $key eq 'compseq'
            ||  $key eq 'seqref'
            ||  $key eq 'featref'
                ) {
           my $clean_id = getCleanID($value);
           $value = $clean_id;
        }
        
        $self->{ 'attr' }->{ $key } = $value;
    }
  }

=item $elem->setattr( $key, $value )

B<Description:> sets an attribute (equivalent to addattr(), but kept for consistency with other API classes)

B<Parameters:>  ($key, $value) - a key value pair 

B<Returns:> None

=cut

sub setattr
  {
      $logger->debug("") if $logger->is_debug;
    my $self = shift;
    my ( $key, $value ) = @_;

    $self->addattr( $key, $value );
  }

=item $elem->setattrh( $hashref )

B<Description:> sets an elements attributes according to the key/value pairs defined in the hash pointed to by $hashref

B<Parameters:> ($hashref) a reference pointing to a hash containing the key value pairs to be added as element attributes

B<Returns:> None

=cut

sub setattrh
  {
      $logger->debug("") if $logger->is_debug;
    my $self = shift;

    my ($href) = @_;

    foreach my $key (keys(%{$href}))
      {
    $self->addattr( $key, $href->{$key} );
      }
  }

=item $elem->dropattr( $key )

B<Description:> removes an attribute from a BSML element

B<Parameters:> ($key) - the attribute key

B<Returns:> None

=cut

sub dropattr
{
      $logger->debug("") if $logger->is_debug;
 
    my $self = shift;
    my ( $key, $value ) = @_;

    delete($self->{'attr'}->{$key});
  }

=item $elem->returnattrHashR()

B<Description:> returns a hash reference to the key, value pairs making up an elements attribute set

B<Parameters:> None

B<Returns:> hash reference

=cut

sub returnattrHashR
  {    
       $logger->debug("") if $logger->is_debug;

   my $self = shift;
    return $self->{'attr'};
  }

=item $elem->returnattr( $key )

B<Description:> returns an attributes value given its key

B<Parameters:> ( $key ) - an attribute key

B<Returns:> an attribute value

=cut

sub returnattr
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ($key) = @_;

    return $self->{'attr'}->{$key};
  }

=item $elem->addBsmlAttr( $key, $value )

B<Description:> add a BSML attribute to a Bsml element

B<Parameters:> ($key, $value) - a key value pair 

B<Returns:> None

=cut

sub addBsmlAttr
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ( $key, $value ) = @_;
    if( defined( $value ) ){
    push(@{$self->{ 'BsmlAttr' }->{ $key }}, $value);
    }
  }

=item $elem->setattr( $key, $value )

B<Description:> sets a BSML attribute (equivalent to addBsmlAttr(), but kept for consistency with other API classes)

B<Parameters:>  ($key, $value) - a key value pair 

B<Returns:> None

=cut

sub setBsmlAttr
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ( $key, $value ) = @_;
    
    $self->addBsmlAttr( $key, $value );
  }

=item $elem->setBsmlAttrh( $hashref )

B<Description:> sets BSML attributes according to the key/value pairs defined in the hash pointed to by $hashref

B<Parameters:> ($hashref) a reference pointing to a hash containing the key value pairs to be added as Bsml attributes

B<Returns:> None

sub setBsmlAttrh
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;

    my ($href) = @_;

    foreach my $key (keys(%{$href}))
      {
    $self->addBsmlAttr( $key, $href->{$key} );
      }
  }


=item $elem->dropBsmlAttr( $key )

B<Description:> removes a BSML attribute from a BSML element

B<Parameters:> ($key) - the attribute key

B<Returns:> None

=cut

sub dropBsmlAttr
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ( $key, $value ) = @_;

    delete($self->{ 'BsmlAttr' }->{ $key });
  }

=item $elem->returnBsmlAttrHashR()

B<Description:> returns a hash reference to the key, value pairs making up an elements attribute set

B<Parameters:> None

B<Returns:> hash reference

=cut

sub returnBsmlAttrHashR
  {
      $logger->debug("") if $logger->is_debug;


    my $self = shift;
    return $self->{ 'BsmlAttr' };
  }

=item $elem->returnBsmlAttr( $key )

B<Description:> returns a BSML attribute value given its key

B<Parameters:> ($key) - the attribute key

B<Returns:> attribute value

=cut

sub returnBsmlAttr
{
      $logger->debug("") if $logger->is_debug;

  my $self = shift;
  my ($key) = @_;
  return $self->{'BsmlAttr'}->{$key};
}

sub addBsmlLink
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ($rel, $href, $role) = @_;

    if( defined($rel) && defined($href) )
    {
        if (defined $role) {
            push(@{$self->{'BsmlLink'}}, {rel=>$rel, href=>$href, role=>$role});
        } else {
            push(@{$self->{'BsmlLink'}}, {rel=>$rel, href=>$href});
        }
    return @{$self->{'BsmlLink'}} - 1;
    }
  }
  
sub setBsmlLink
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ($rel, $href) = @_;
    
    return $self->addBsmlLink( $rel, $href );
  }

sub dropBsmlLink
  {
      $logger->debug("") if $logger->is_debug;


    my $self = shift;
    my ($index) = @_;

    my $newlist;

    for( my $i=0; $i< @{$self->{'BsmlLink'}}; $i++ )
      {
    if( $i != $index )
      {
        push( @{$newlist}, $self->{'BsmlLink'}[$i] );
      }
      }

    $self->{'BsmlLink'} = $newlist;
  }

sub returnBsmlLinkListR
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    return $self->{'BsmlLink'};
  }

sub returnBsmlLinkR
  {
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ($index) = @_;

    return $self->{'BsmlLink'}[$index];
  }

#----------------------------------------------------------------------
# BsmlCrossReference support
#
#----------------------------------------------------------------------

sub addBsmlCrossReference
{
      $logger->debug("") if $logger->is_debug;


    my $self = shift;

    push (@{$self->{'BsmlCrossReference'}}, new BSML::BsmlCrossReference);
    my $index = @{$self->{'BsmlCrossReference'}} - 1;

    #print Dumper $self->{'BsmlCrossReference'}[$index];
    #print $index;
    return $index;
}

sub returnBsmlCrossReferenceR
{
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ($index) = @_;

    return $self->{'BsmlCrossReference'}[$index];
}

sub returnBsmlCrossReferenceListR {

      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    return $self->{'BsmlCrossReference'};
    
}

sub dropBsmlCrossReference
{
      $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my $index = shift;

    delete($self->{'BsmlCrossReference'}->[$index]);
}

sub write()
  {
      $logger->debug("") if $logger->is_debug;
    
  }


#
# 2004-11-04 Support to add, drop, return BsmlAttributeList 
# 

sub addBsmlAttributeList {


    $logger->debug("") if $logger->is_debug;

    my $self = shift;
    my ( $listref ) = @_;

    push (@{$self->{'BsmlAttributeList'}}, $listref);
    
}

sub dropBsmlAttributeList  {

      $logger->debug("") if $logger->is_debug;

    my $self = shift;

    $self->{ 'BsmlAttributeList' } = undef;
  }


sub returnBsmlAttributeListHashR {

    $logger->debug("") if $logger->is_debug;
    
    my $self = shift;
    
#    $logger->fatal(Dumper $self);

    return $self->{ 'BsmlAttributeList' };

}

## getCleanID is used to convert strings into XML 1.0 compliant id strings
## see: http://www.w3.org/TR/REC-xml/#NT-Letter
sub getCleanID {
    my ($self, $text);
    
    my $basechar =  '\x{0041}-\x{005A}\x{0061}-\x{007A}\x{00C0}-\x{00D6}\x{00D8}-\x{00F6}'
                   .'\x{00F8}-\x{00FF}\x{0100}-\x{0131}\x{0134}-\x{013E}\x{0141}-\x{0148}'
                   .'\x{014A}-\x{017E}\x{0180}-\x{01C3}\x{01CD}-\x{01F0}\x{01F4}-\x{01F5}'
                   .'\x{01FA}-\x{0217}\x{0250}-\x{02A8}\x{02BB}-\x{02C1}\x{0386}\x{0388}-\x{038A}'
                   .'\x{038C}\x{038E}-\x{03A1}\x{03A3}-\x{03CE}\x{03D0}-\x{03D6}\x{03DA}'
                   .'\x{03DC}\x{03DE}\x{03E0}\x{03E2}-\x{03F3}\x{0401}-\x{040C}\x{040E}-\x{044F}'
                   .'\x{0451}-\x{045C}\x{045E}-\x{0481}\x{0490}-\x{04C4}\x{04C7}-\x{04C8}'
                   .'\x{04CB}-\x{04CC}\x{04D0}-\x{04EB}\x{04EE}-\x{04F5}\x{04F8}-\x{04F9}'
                   .'\x{0531}-\x{0556}\x{0559}\x{0561}-\x{0586}\x{05D0}-\x{05EA}\x{05F0}-\x{05F2}'
                   .'\x{0621}-\x{063A}\x{0641}-\x{064A}\x{0671}-\x{06B7}\x{06BA}-\x{06BE}'
                   .'\x{06C0}-\x{06CE}\x{06D0}-\x{06D3}\x{06D5}\x{06E5}-\x{06E6}\x{0905}-\x{0939}'
                   .'\x{093D}\x{0958}-\x{0961}\x{0985}-\x{098C}\x{098F}-\x{0990}\x{0993}-\x{09A8}'
                   .'\x{09AA}-\x{09B0}\x{09B2}\x{09B6}-\x{09B9}\x{09DC}-\x{09DD}\x{09DF}-\x{09E1}'
                   .'\x{09F0}-\x{09F1}\x{0A05}-\x{0A0A}\x{0A0F}-\x{0A10}\x{0A13}-\x{0A28}'
                   .'\x{0A2A}-\x{0A30}\x{0A32}-\x{0A33}\x{0A35}-\x{0A36}\x{0A38}-\x{0A39}'
                   .'\x{0A59}-\x{0A5C}\x{0A5E}\x{0A72}-\x{0A74}\x{0A85}-\x{0A8B}\x{0A8D}'
                   .'\x{0A8F}-\x{0A91}\x{0A93}-\x{0AA8}\x{0AAA}-\x{0AB0}\x{0AB2}-\x{0AB3}'
                   .'\x{0AB5}-\x{0AB9}\x{0ABD}\x{0AE0}\x{0B05}-\x{0B0C}\x{0B0F}-\x{0B10}'
                   .'\x{0B13}-\x{0B28}\x{0B2A}-\x{0B30}\x{0B32}-\x{0B33}\x{0B36}-\x{0B39}'
                   .'\x{0B3D}\x{0B5C}-\x{0B5D}\x{0B5F}-\x{0B61}\x{0B85}-\x{0B8A}\x{0B8E}-\x{0B90}'
                   .'\x{0B92}-\x{0B95}\x{0B99}-\x{0B9A}\x{0B9C}\x{0B9E}-\x{0B9F}\x{0BA3}-\x{0BA4}'
                   .'\x{0BA8}-\x{0BAA}\x{0BAE}-\x{0BB5}\x{0BB7}-\x{0BB9}\x{0C05}-\x{0C0C}'
                   .'\x{0C0E}-\x{0C10}\x{0C12}-\x{0C28}\x{0C2A}-\x{0C33}\x{0C35}-\x{0C39}'
                   .'\x{0C60}-\x{0C61}\x{0C85}-\x{0C8C}\x{0C8E}-\x{0C90}\x{0C92}-\x{0CA8}'
                   .'\x{0CAA}-\x{0CB3}\x{0CB5}-\x{0CB9}\x{0CDE}\x{0CE0}-\x{0CE1}\x{0D05}-\x{0D0C}'
                   .'\x{0D0E}-\x{0D10}\x{0D12}-\x{0D28}\x{0D2A}-\x{0D39}\x{0D60}-\x{0D61}'
                   .'\x{0E01}-\x{0E2E}\x{0E30}\x{0E32}-\x{0E33}\x{0E40}-\x{0E45}\x{0E81}-\x{0E82}'
                   .'\x{0E84}\x{0E87}-\x{0E88}\x{0E8A}\x{0E8D}\x{0E94}-\x{0E97}\x{0E99}-\x{0E9F}'
                   .'\x{0EA1}-\x{0EA3}\x{0EA5}\x{0EA7}\x{0EAA}-\x{0EAB}\x{0EAD}-\x{0EAE}'
                   .'\x{0EB0}\x{0EB2}-\x{0EB3}\x{0EBD}\x{0EC0}-\x{0EC4}\x{0F40}-\x{0F47}'
                   .'\x{0F49}-\x{0F69}\x{10A0}-\x{10C5}\x{10D0}-\x{10F6}\x{1100}\x{1102}-\x{1103}'
                   .'\x{1105}-\x{1107}\x{1109}\x{110B}-\x{110C}\x{110E}-\x{1112}\x{113C}'
                   .'\x{113E}\x{1140}\x{114C}\x{114E}\x{1150}\x{1154}-\x{1155}\x{1159}'
                   .'\x{115F}-\x{1161}\x{1163}\x{1165}\x{1167}\x{1169}\x{116D}-\x{116E}'
                   .'\x{1172}-\x{1173}\x{1175}\x{119E}\x{11A8}\x{11AB}\x{11AE}-\x{11AF}'
                   .'\x{11B7}-\x{11B8}\x{11BA}\x{11BC}-\x{11C2}\x{11EB}\x{11F0}\x{11F9}'
                   .'\x{1E00}-\x{1E9B}\x{1EA0}-\x{1EF9}\x{1F00}-\x{1F15}\x{1F18}-\x{1F1D}'
                   .'\x{1F20}-\x{1F45}\x{1F48}-\x{1F4D}\x{1F50}-\x{1F57}\x{1F59}\x{1F5B}'
                   .'\x{1F5D}\x{1F5F}-\x{1F7D}\x{1F80}-\x{1FB4}\x{1FB6}-\x{1FBC}\x{1FBE}'
                   .'\x{1FC2}-\x{1FC4}\x{1FC6}-\x{1FCC}\x{1FD0}-\x{1FD3}\x{1FD6}-\x{1FDB}'
                   .'\x{1FE0}-\x{1FEC}\x{1FF2}-\x{1FF4}\x{1FF6}-\x{1FFC}\x{2126}\x{212A}-\x{212B}'
                   .'\x{212E}\x{2180}-\x{2182}\x{3041}-\x{3094}\x{30A1}-\x{30FA}\x{3105}-\x{312C}'
                   .'\x{AC00}-\x{D7A3}';
    my $ideographic = '\x{4E00}-\x{9FA5}\x{3007}\x{3021}-\x{3029}';
    my $letter = $basechar . $ideographic;
    my $digit =  '\x{0030}-\x{0039}\x{0660}-\x{0669}\x{06F0}-\x{06F9}\x{0966}-\x{096F}'
                .'\x{09E6}-\x{09EF}\x{0A66}-\x{0A6F}\x{0AE6}-\x{0AEF}\x{0B66}-\x{0B6F}'
                .'\x{0BE7}-\x{0BEF}\x{0C66}-\x{0C6F}\x{0CE6}-\x{0CEF}\x{0D66}-\x{0D6F}'
                .'\x{0E50}-\x{0E59}\x{0ED0}-\x{0ED9}\x{0F20}-\x{0F29}';
    my $combiningchar = '\x{0300}-\x{0345}\x{0360}-\x{0361}\x{0483}-\x{0486}\x{0591}-\x{05A1}'
                .'\x{05A3}-\x{05B9}\x{05BB}-\x{05BD}\x{05BF}\x{05C1}-\x{05C2}\x{05C4}'
                .'\x{064B}-\x{0652}\x{0670}\x{06D6}-\x{06DC}\x{06DD}-\x{06DF}\x{06E0}-\x{06E4}'
                .'\x{06E7}-\x{06E8}\x{06EA}-\x{06ED}\x{0901}-\x{0903}\x{093C}\x{093E}-\x{094C}'
                .'\x{094D}\x{0951}-\x{0954}\x{0962}-\x{0963}\x{0981}-\x{0983}\x{09BC}\x{09BE}'
                .'\x{09BF}\x{09C0}-\x{09C4}\x{09C7}-\x{09C8}\x{09CB}-\x{09CD}\x{09D7}'
                .'\x{09E2}-\x{09E3}\x{0A02}\x{0A3C}\x{0A3E}\x{0A3F}\x{0A40}-\x{0A42}'
                .'\x{0A47}-\x{0A48}\x{0A4B}-\x{0A4D}\x{0A70}-\x{0A71}\x{0A81}-\x{0A83}\x{0ABC}'
                .'\x{0ABE}-\x{0AC5}\x{0AC7}-\x{0AC9}\x{0ACB}-\x{0ACD}\x{0B01}-\x{0B03}\x{0B3C}'
                .'\x{0B3E}-\x{0B43}\x{0B47}-\x{0B48}\x{0B4B}-\x{0B4D}\x{0B56}-\x{0B57}'
                .'\x{0B82}-\x{0B83}\x{0BBE}-\x{0BC2}\x{0BC6}-\x{0BC8}\x{0BCA}-\x{0BCD}'
                .'\x{0BD7}\x{0C01}-\x{0C03}\x{0C3E}-\x{0C44}\x{0C46}-\x{0C48}\x{0C4A}-\x{0C4D}'
                .'\x{0C55}-\x{0C56}\x{0C82}-\x{0C83}\x{0CBE}-\x{0CC4}\x{0CC6}-\x{0CC8}'
                .'\x{0CCA}-\x{0CCD}\x{0CD5}-\x{0CD6}\x{0D02}-\x{0D03}\x{0D3E}-\x{0D43}'
                .'\x{0D46}-\x{0D48}\x{0D4A}-\x{0D4D}\x{0D57}\x{0E31}\x{0E34}-\x{0E3A}'
                .'\x{0E47}-\x{0E4E}\x{0EB1}\x{0EB4}-\x{0EB9}\x{0EBB}-\x{0EBC}\x{0EC8}-\x{0ECD}'
                .'\x{0F18}-\x{0F19}\x{0F35}\x{0F37}\x{0F39}\x{0F3E}\x{0F3F}\x{0F71}-\x{0F84}'
                .'\x{0F86}-\x{0F8B}\x{0F90}-\x{0F95}\x{0F97}\x{0F99}-\x{0FAD}\x{0FB1}-\x{0FB7}'
                .'\x{0FB9}\x{20D0}-\x{20DC}\x{20E1}\x{302A}-\x{302F}\x{3099}\x{309A}';
    my $extender = '\x{00B7}\x{02D0}\x{02D1}\x{0387}\x{0640}\x{0E46}\x{0EC6}\x{3005}'
                  .'\x{3031}-\x{3035}\x{309D}-\x{309E}\x{30FC}-\x{30FE}';
    
    ## $namechar = $letter . $digit . '\.' . '\-' . '_' . ':' . $combiningchar . $extender;
    
    ## support invoking as an object or class method
    if (eval{$_[0]->isa('BSML::BsmlElement')}) {
        ($self, $text)  = @_;
    } else {
        ($text)         = @_;
    }

    ## enforce XML spec on first char   
    if ($text =~ /^[^$letter\_:]{1}/) {
        $text = '_'.$text;
    }
    
    ## enforce XML spec on remaining chars
    $text =~ s/[^$letter$digit\.\-_:$combiningchar$extender]{1}/_/g;

    return $text;
}


1
