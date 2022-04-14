# HGTfinder
Workflow of HGT identification:

Step 1:
AI calculation (see example in example4AI_estimate folder)

perl calculate_alien_index_MD_Xing.pl -c Spathaspora_gorwiae_ai.config

Step 2:
building tree:

/tools/mafft-7.299-without-extensions/bin/mafft --thread 6 --auto $filename > ${gene}_aln.fasta
/tools/trimal-trimAl/source/trimal -in ${gene}_aln.fasta -out ${gene}_aln_trimmed.fasta -automated1
/tools/IQ-TREE-1.5.1/build/iqtree-omp -nt 6 -st AA -s ${gene}_aln_trimmed.fasta -m TEST -mrate G4 -keep-ident -bb 1000 -pre ${gene}

Midpoint root a phylogeny with R:

mytree <- read.tree("$filename")
midponit_mytree <- ladderize (midpoint(mytree))
write.tree (midponit_mytree, "${gene}midpoint.tree")

For questions or suggestions, please contact xingxingshen@zju.edu.cn

When using the pipelibe in published research, please cite:

Shen X.-X., Opulente D.A., Kominek J.#, Zhou X., Steenwyk J., Buh K.V., Haase M., Wisecaver J.H., Wang M., Doering D.T., Boudouris J., Schneider R., Langdon Q.K., Ohkuma M., Endoh R., Takashima M., Manabe R., Čadež N., Libkind D., Rosa C., DeVirgilio J., Hulfachor A., Groenewald M., Kurtzman C.P., Hittinger C.T., Rokas A. 2018. The tempo and mode of genome evolution in the budding yeast subphylum. Cell 175: 1533-1545. doi:10.1016/j.cell.2018.10.023
