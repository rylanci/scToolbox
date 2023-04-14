This repository documents the tools and pipelines our bioinformatic team apply to single cell datasets. 
The goal is for it to contain templates of our analysis that are broadly applicable to data of various tissues with mouse and human support. 
Ideally this will be useful to us for version control and also make these scripts more accessable for anybody with a use for them:
whether that means using them as is or making them your own. 

###### TO DO:

	- make jupyter-nb tutorials in each of the toolsets
	
	- generate test datasets to process through each toolset
		1. a small object (rna and atac)
		2. bulk matrices pertaining to that object  
		3. small fastqs for sci-ATAC
		4. small filter feature matrix for RNA


	- create yamls for all conda env necessary to run scripts in this folder 
		make sure it can be loaded via conda yaml thing. yaml can then be provided


	- homepage markdown could explain the schematic and guide towards jupyter tutorials
		explain how one could create env using yaml & conda. Run on example data.
		then transfer to additional data 
