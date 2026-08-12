use strict;
use warnings;

use Bio::EnsEMBL::Registry;

Bio::EnsEMBL::Registry->load_registry_from_url(
  'mysql://ensro@mysql-ens-havana-prod-1:4581/havana_human_t2t?group=core&species=homo_sapiens_gca009914755v4'
);

1;
