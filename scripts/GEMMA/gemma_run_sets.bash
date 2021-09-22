#!/usr/bin/env bash

PLINK_INPUT="plink.nam_stringent" # The output name for Plink data
DESIGN_FILE="gemma_design.csv"
RUNS=$(wc -l ${DESIGN_FILE} | cut -f 1 -d ' ')
KINSHIP="s"  # c=centered s=standardized

for i in $(seq 1 ${RUNS})
do
    DESIGN=$(cat ${DESIGN_FILE} | head -n ${i} | tail -n 1)
    IFS=',' read -ra ARRAY <<< "${DESIGN}"
    RUN_NAME=${ARRAY[0]}
    GEMMA_OUTPUT=${PLINK_INPUT}_${RUN_NAME}
    Y_VARS=${ARRAY[1]}
    COVARS=${ARRAY[2]}
    echo ${COVARS}

    if [ "${COVARS}" != "0" ]
    then
        ~/gemma-0.98.3 \
            -bfile ${PLINK_INPUT} \
            -k ./output/${PLINK_INPUT}.${KINSHIP}XX.txt \
            -lmm 1 \
            -o ${GEMMA_OUTPUT} \
            -n ${Y_VARS} \
            -c ${COVARS}
    fi

    if [ "${COVARS}" == "0" ]
    then
        ~/gemma-0.98.3 \
            -bfile ${PLINK_INPUT} \
            -k ./output/${PLINK_INPUT}.${KINSHIP}XX.txt \
            -lmm 1 \
            -o ${GEMMA_OUTPUT} \
            -n ${Y_VARS} 
    fi
done

# FAM 6, 7, 8, 9, 10, 11, 12, 13 ### THIS IS WRONG>>>
# maturity, height, biomass, DTH, MeterWeight, StandCount, DryWeight, WetWeight
# FAM 13-46
# ADF,AD.ICP,Adj_CP,aNDF,aNDFom,Ash,Ca,Cl,Crude.protein,DCAD,Dry.Matter,EE.Fat,K,Lignin,Lignin_.NDF,Mg,Moisture,Na,NDF,NDICP_w.oNa2SO3,NEG_OARDC,NEL3x_ADF,NEL3x_OARDC,NEM_ADF,NEM_OARDC,NFC,P,RFV,S,SP.CP,Starch,TDN_ADF,TDN_OARDC,WSC_Sugar
# 1 5 7 8
# 1 2 3 4
# ~/gemma-0.98.3 -bfile ${PLINK_OUTPUT} \
#     -k ./output/${PLINK_OUTPUT}.sXX.txt \
#     -lmm 4 -n 1 2 3 4 -o "${PLINK_OUTPUT}_MV"
