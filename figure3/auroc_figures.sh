# Fully normalized cross-platform
python cross_platform_auroc.py --pl ILLUMINA -p AF -c platforms -k line -n --normalize-tp $COMBINED_RESULTS plots/by_platform_fully_normalized/pl_af_illumina_all_normalized_roc.png

python cross_platform_auroc.py --pl OXFORD_NANOPORE -p AF -c platforms -k line -n --normalize-tp $COMBINED_RESULTS plots/by_platform_fully_normalized/pl_af_ont_all_normalized_roc.png

python cross_platform_auroc.py --pl ILLUMINA -p DP_ALT -c platforms -k line -n --normalize-tp $COMBINED_RESULTS plots/by_platform_fully_normalized/pl_dpalt_illumina_all_normalized_roc.png

python cross_platform_auroc.py --pl OXFORD_NANOPORE -p DP_ALT -c platforms -k line -n --normalize-tp $COMBINED_RESULTS plots/by_platform_fully_normalized/pl_dpalt_ont_all_normalized_roc.png

# Fully normalized cross-groups
python cross_platform_auroc.py --pl ILLUMINA -p AF -c groups -k line -n --normalize-tp $COMBINED_RESULTS plots/by_group_fully_normalized/gr_af_illumina_all_normalized_roc.png

python cross_platform_auroc.py --pl OXFORD_NANOPORE -p AF -c groups -k line -n --normalize-tp $COMBINED_RESULTS plots/by_group_fully_normalized/gr_af_ont_all_normalized_roc.png

python cross_platform_auroc.py --pl ILLUMINA -p DP_ALT -c groups -k line -n --normalize-tp $COMBINED_RESULTS plots/by_group_fully_normalized/gr_dpalt_illumina_all_normalized_roc.png

python cross_platform_auroc.py --pl OXFORD_NANOPORE -p DP_ALT -c groups -k line -n --normalize-tp $COMBINED_RESULTS plots/by_group_fully_normalized/gr_dpalt_ont_all_normalized_roc.png
