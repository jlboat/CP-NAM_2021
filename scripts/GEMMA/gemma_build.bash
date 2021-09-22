#!/usr/bin/env bash

VCF="NAM.imputed.ann.vcf"
PLINK_OUTPUT="plink.nam_stringent" # The output name for Plink data
#PHENOTYPES="mean_ama_BLUP_tan.txt" # tab-separate file with header = FID\tIID\tPhenotype1\tPhenotype2\t...PhenotypeN
PHENOTYPES="CP-NAM_NIR_MV.gemma.tsv"
MISSING=0.3
HWE=0
MAF=0.05
STANDARDIZED_KINSHIP="TRUE"
EIGEN="TRUE"

~/Plink/plink \
  --double-id \
  --allow-no-sex \
  --all-pheno \
  --make-bed \
  --out ${PLINK_OUTPUT} \
  --pheno ${PHENOTYPES} \
  --vcf ${VCF}

if [ ${STANDARDIZED_KINSHIP} == "TRUE" ]
then
    # gk=1 (centered relatedness matrix) gk=2 (standardized relatedness matrix)
    /zfs/tillers/panicle/software/gemma-0.98.3 -bfile ${PLINK_OUTPUT} -gk 2 \
       -o ${PLINK_OUTPUT} -miss ${MISSING} -hwe ${HWE} -maf ${MAF}
else
    # gk=1 (centered relatedness matrix) gk=2 (standardized relatedness matrix)
    /zfs/tillers/panicle/software/gemma-0.98.3 -bfile ${PLINK_OUTPUT} -gk 1 \
       -o ${PLINK_OUTPUT} -miss ${MISSING} -hwe ${HWE} -maf ${MAF}
fi

if [ ${EIGEN} == "TRUE" ]
then
    if [ ${STANDARDIZED_KINSHIP} == "TRUE" ]
    then
        # Run gemma eigenvalue decomposition
        /zfs/tillers/panicle/software/gemma-0.98.3 -bfile ${PLINK_OUTPUT} \
            -k ./output/${PLINK_OUTPUT}.sXX.txt \
            -eigen -o ${PLINK_OUTPUT} \
            -miss ${MISSING} -hwe ${HWE} -maf ${MAF}
    else
        # Run gemma eigenvalue decomposition
        /zfs/tillers/panicle/software/gemma-0.98.3 -bfile ${PLINK_OUTPUT} \
            -k ./output/${PLINK_OUTPUT}.cXX.txt \
            -eigen -o ${PLINK_OUTPUT} \
            -miss ${MISSING} -hwe ${HWE} -maf ${MAF}
    fi
fi
