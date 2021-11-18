###### script to set jobs going
# submit for each chromosome


gitdir=/net/harris/vol1/home/beichman/scriptsAndGitDirs/DNA_Shape_Project
scriptdir=$gitdir/analyses/mutyper/apes/
script=$scriptdir/runMutyper.AnyVCFREF.NOSTRICT.apes.REDOKSFS.sh
# originally in /net/harris/vol1/project/primate_ervs/hg18_references_with_gagp_ancestral_alleles_exclude_recurrent
# but I copied to : so I can make .fai files 
# cp /net/harris/vol1/project/primate_ervs/hg18_references_with_gagp_ancestral_alleles_exclude_recurrent/*fa 
refdir=/net/harris/vol1/home/beichman/apes/polarized_ref_fastas/hg18_references_with_gagp_ancestral_alleles_exclude_recurrent 


vcfdir=/net/harris/vol1/data/great_ape_genome_project/eichlerlab.gs.washington.edu/greatape/data/VCFs/SNPs

#chimp low cov sample to exclude: Pan_troglodytes_ellioti-Banyo
#gorilla low cov samples to exclude: Gorilla_gorilla_gorilla-X00108_Abe , Gorilla_gorilla_gorilla-KB7973_Porta 


speciesList='Gorilla Pan_paniscus Pan_troglodytes Pongo_abelii Pongo_pygmaeus'
#speciesList="Gorilla"
### need to update name of each fasta file to be chr1, 2 etc.:
# once (only in my dir)
#for chr in {1..22}
#do
#sed -i"" "s/^>ANCESTOR.*/>chr$chr/g" $refdir/homo_sapiens_ancestor_${chr}.fa 
#done

### need to separate out by population
#for chr in {1..22}
#for chr in 1
#for chr in {2..3}
#for chr in 1 2
#for chr in {13..22}

for species in $speciesList
do
echo $species
# set inds to remove
# use = with [] and == with [[]]
if [ $species = "Gorilla" ]
then
indsToRemove='Gorilla_gorilla_gorilla-X00108_Abe,Gorilla_gorilla_gorilla-KB7973_Porta'
elif [ $species = "Pan_troglodytes" ]
then
indsToRemove='Pan_troglodytes_ellioti-Banyo'
else
indsToRemove="none"
fi

for chr in {1..22}
do 

label=chr${chr}
reference=${species}_${label}.fa
vcf=${species}.vcf.gz


label="chr${chr}"
qsub -N ${species}.chr.${chr} $script $vcfdir/$vcf $refdir/$reference $species $label $indsToRemove



done
done
