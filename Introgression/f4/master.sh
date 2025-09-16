P1_list=(isocomma_simulator isocomma_simulator mothone_mothone menophilus_hicetas menophilus_ssp)
P2_list=(isocomma_isocomma isocomma_isocomma mothone_messenina menophilus_ernestoi menophilus_ernestoi)
P3_list=(menophilus_zaneka menophilus_menophilus menophilus_zaneka cydon mothone_mothone)
P4_list=(menophilus_ernestoi menophilus_ernestoi menophilus_ernestoi maeolus mothone_messenina)

for i in ${!P1_list[@]}; do
    P1=${P1_list[$i]}
    P2=${P2_list[$i]}
    P3=${P3_list[$i]}
    P4=${P4_list[$i]}
    
    sbatch ./f4.sh ${P1} ${P2} ${P3} ${P4}
done
