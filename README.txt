Codes for analyzing biological shape using inconsistent surface registration

Reference:
G. P. T. Choi, D. Qiu and L. M. Lui, 
"Shape analysis via inconsistent surface registration."
Proceedings of the Royal Society A, 476(2242), 20200147, 2020.

Copyright (c) 2020, Gary Pui-Tung Choi, Di Qiu, Lok Ming Lui

===============================================================
Usage:
run_pairwise_map.m: Compute a mapping between two molars using our proposed method with comparison
run_pairwise_map_comparison_tmap.m: Compare the performance of surface mapping between our method and the Teichmuller mapping method
run_inconsistent_shape_registration_all.m: Run inconisistent shape registration for the 50x50 pairs of meshes
run_clustering_binary: Cluster the specimens into 2 groups based on our inconsistent surface mapping result
run_clustering_genera: Cluster the specimens into 5 groups based on our inconsistent surface mapping result
run_clustering_genera_comparison: Classification into five groups using traditional methods

The tooth dataset is from https://gaotingran.com/resources/auto3dgm_hdm/