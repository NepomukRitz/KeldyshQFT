#!/bin/sh

# sync Keldysh_mfRG_data between the users Sa.Aguirre, Fabian.Kugler, E.Walter

# files are stored in /project/th-sratch/first-character-of-user-name/user-name/Keldysh_mfRG_data
# user-name is stored as $USER
initusertmp="$(echo $USER | head -c 1)" # obtain the first character
initusertmp="$(echo $initusertmp | tr '[:upper:]' '[:lower:]')" # to convert into lowercase

# make folder if not existing
mkdir -p /project/th-scratch/"$initusertmp"/"$USER"/Keldysh_mfRG_data/

echo "\t\t --- \t Sync Keldysh_mfRG_data \t --- \t\t"

# sync folders from co-workers with own folder
rsync -avuh /project/th-scratch/s/Sa.Aguirre/Keldysh_mfRG_data/ /project/th-scratch/"$initusertmp"/"$USER"/Keldysh_mfRG_data/
rsync -avuh /project/th-scratch/f/Fabian.Kugler/Keldysh_mfRG_data/ /project/th-scratch/"$initusertmp"/"$USER"/Keldysh_mfRG_data/
rsync -avuh /project/th-scratch/e/E.Walter/Keldysh_mfRG_data/ /project/th-scratch/"$initusertmp"/"$USER"/Keldysh_mfRG_data/

echo "\t\t --- \t Keldysh_mfRG_data synced \t --- \t\t"
