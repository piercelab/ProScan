#!/bin/bash
PROJECT_DIR=$(dirname "${BASH_SOURCE[0]}")

####### SETTINGS #######
# Location of your bash script for cron to run:
UPDATE_SCRIPT=$PROJECT_DIR/updates/cronscript.sh
# Set NEED_CONDA to 0 if the update script involves python scripts.
NEED_CONDA=true


####### LOGIC #######
# Quit without error if script doesn't exist.
[[ -e $UPDATE_SCRIPT ]] || exit

SCR_DIR=$(dirname "$UPDATE_SCRIPT")
SCR_NAME=$(basename "$UPDATE_SCRIPT")

# Move into script directory
cd "$SCR_DIR"
# Run with conda if needed, or with bash if not.
if $NEED_CONDA ; then
    # Source conda definitions in case cron is missing the commands.
    source "/etc/profile.d/conda.sh"
    # Run the script in the conda environment
    conda run -n wsgi bash "$SCR_NAME"
else
    bash "$SCR_NAME"
fi