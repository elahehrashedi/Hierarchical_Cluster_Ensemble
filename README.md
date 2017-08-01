# A Hierarchical Clusterer Ensemble Method Based on Boosting Theory

## Citation

Please cite this work in your publications if it helps your research:

	@article{rashedi2013hierarchical,
		title={A hierarchical clusterer ensemble method based on boosting theory},
		author={Rashedi, Elaheh and Mirzaei, Abdolreza},
		journal={Knowledge-Based Systems},
		volume={45},
		pages={83--93},
		year={2013},
		publisher={Elsevier}
	}

## Abstract:

Bagging and boosting are two well-known methods of developing classifier ensembles. It is generally agreed that the clusterer ensemble methods that utilize the boosting concept can create clusterings with quality and robustness improvements. In this paper, we introduce a new boosting based hierarchical clusterer ensemble method called Bob-Hic. This method is utilized to create a consensus hierarchical clustering (h-clustering) on a dataset, which is helpful to improve the clustering accuracy. Bob-Hic includes several boosting iterations. In each iteration, first, a weighted random sampling is performed on the original dataset. An individual h-clustering is then created on the selected samples. At the end of the iterations, the individual clusterings are combined to a final consensus h-clustering. The intermediate structures used in the combination are distance descriptor matrices which correspond to individual h-clustering results. This final integration is done through an information theoretic approach. Experiments on popular synthetic and real datasets confirm that the proposed method improves the results of simple clustering algorithms. In addition, our experimental results confirm that this method provides better consensus clustering quality compared to other available ensemble techniques.