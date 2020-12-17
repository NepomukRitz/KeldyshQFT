#!/bin/sh

# sync Keldysh_mfRG_data between the users Sa.Aguirre, Fabian.Kugler, E.Walter

# files are stored in /project/th-sratch/first-character-of-user-name/user-name/Keldysh_mfRG_data
# user-name is stored as $USER
initusertmp="$(echo $USER | head -c 1)" # obtain the first character
initusertmp="$(echo $initusertmp | tr '[:upper:]' '[:lower:]')" # to convert into lowercase

# make exchange folder if not existing
mkdir -p /project/th-scratch/"$initusertmp"/"$USER"/Keldysh_mfRG_data/

# set rights for exchange folder (recursively) using 750: full rights for user, read and execute rights for group, no rights for others
chmod -R 750 /project/th-scratch/"$initusertmp"/"$USER"/Keldysh_mfRG_data/

echo "\t\t --- \t Start sync of Keldysh_mfRG_data \t --- \t\t"

# sync folders from co-workers with own folder
echo "\t\t --- \t Sync data from Sa.Aguirre \t --- \t\t"
rsync -avuh /project/th-scratch/s/Sa.Aguirre/Keldysh_mfRG_data/ /project/th-scratch/"$initusertmp"/"$USER"/Keldysh_mfRG_data/
echo "\t\t --- \t Sync data from Fabian.Kugler \t --- \t\t"
rsync -avuh /project/th-scratch/f/Fabian.Kugler/Keldysh_mfRG_data/ /project/th-scratch/"$initusertmp"/"$USER"/Keldysh_mfRG_data/
echo "\t\t --- \t Sync data from E.Walter \t --- \t\t"
rsync -avuh /project/th-scratch/e/E.Walter/Keldysh_mfRG_data/ /project/th-scratch/"$initusertmp"/"$USER"/Keldysh_mfRG_data/

echo "\t\t --- \t Finish sync of Keldysh_mfRG_data \t --- \t\t"
