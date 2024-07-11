This 'desi_release' branch is used to generate DLA catalogs for DESI releases.

Before running the GP DLA finder, you need to prepare preprocess files and trained model.
For preprocess files, you can find how to generate it from preprocessing procedures. (Link:)
All other training data are provided in /data directory except the catalog which is too large to upload. You can download the catalog file here: https://drive.google.com/drive/folders/1gPm-pJPbc_SpWdI_6qnUrRJUN5PbUlUu?usp=sharing.

Run script_for_gp.m in code directory to process each healpix file.

Run generate_catalog.py to generate DLA catalogs for each healpix and then combine them together.
