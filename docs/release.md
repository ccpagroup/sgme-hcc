To create a release zip
-----------------------
```sh
# Get today's date
today=$(date +%y%m%d-%H%M)

# Create the archive file with today's date as the file name
git archive --prefix sgme-hcc/ --format=tar main | bzip2 > release/sgme_$today.bz2
```