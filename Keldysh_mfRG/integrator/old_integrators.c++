#include "old_integrators.hpp"

FrequencyGrid<eliasGrid> frequencyGrid_bos ('b', 1, Lambda_ini);
FrequencyGrid<eliasGrid> frequencyGrid_fer ('f', 1, Lambda_ini);
rvec bfreqs = frequencyGrid_bos.get_all_frequencies();
rvec ffreqs = frequencyGrid_fer.get_all_frequencies();

