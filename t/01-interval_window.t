use strict;
use warnings;

use Test::More;
use AlignDB::IntSpan;

use AlignDB::Window;

{
    print "#interval_window\n";

    my $maker = AlignDB::Window->new;

    my @data = (
        [   [ AlignDB::IntSpan->new->add_pair( 1, 1000 ), 1, 1000 ],
            [   {   snp_all_bases   => "GT",
                    snp_mutant_to   => "G<->T",
                    snp_query_base  => "T",
                    snp_freq        => 1,
                    snp_occured     => 10,
                    snp_pos         => 9,
                    snp_target_base => "G",
                },
            ],
        ],
    );

    for my $i ( 0 .. $#data ) {
        my ( $input_refs, $except_ref ) = @{ $data[$i] };

        my @results = $maker->interval_window( @{$input_refs} );
        is_deeply( \@results, $except_ref, "interval window $i" );
    }
}
