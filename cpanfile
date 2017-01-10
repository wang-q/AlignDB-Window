requires 'Moose';
requires 'List::Util';
requires 'YAML::Syck';
requires 'AlignDB::IntSpan';
requires 'perl', '5.008001';

on test => sub {
    requires 'Test::More', 0.88;
};
