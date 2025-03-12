#!/bin/bash

# Remove Zone.Identifier files from the current directory and all subdirectories
# These files are created by Windows when downloading files from the internet
find . -type f -name "*:Zone.Identifier" -exec rm {} \;
