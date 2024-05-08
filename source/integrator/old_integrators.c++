#include "old_integrators.hpp"

FrequencyGrid<eliasGrid> frequencyGrid_bos ('b', 1, Lambda_ini, fRG_config());
FrequencyGrid<eliasGrid> frequencyGrid_fer ('f', 1, Lambda_ini, fRG_config());
vec<freqType> bfreqs = frequencyGrid_bos.get_all_frequencies();
vec<freqType> ffreqs = frequencyGrid_fer.get_all_frequencies();

