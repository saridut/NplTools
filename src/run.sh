#!/usr/bin/env bash

dir="/Volumes/lasdata/npl_8x8/mch"

fn_cfg="${dir}/slc_ML4_8x8_45.lmp" #Configuration file name
fn_molinfo="${dir}/slc_ML4_8x8_45_mg.txt" #Molinfo file name

fn_fns=${dir}/traj/markers #Name of file containing dump file names & markers
#fn_fns='dmp_fns' #Name of file containing dump file names

#fn_dump="${dir}/trajeq_lig.bin" #Dump file name/pattern

#Write dump file names to `fn_fns`
#if [[ "${fn_dump}" != *'*'* ]]; then
#    if [[ -f "${fn_dump}" ]]; then 
#        fns=(${fn_dump})
#    else
#        printf "%s\n" "Dump file(s) do not exist. Exiting ..."
#        exit 1
#    fi
#else
#    IFS='*'; read -ra tmp <<< "$fn_dump"; unset IFS
#    prefix="${tmp[0]}"; suffix="${tmp[1]}"
#    np=${#prefix}; ns=${#suffix}
#    fns=()
#    for each in "${prefix}"*"${suffix}"; do
#        if [[ -f "${each}" ]]; then fns+=("${each}"); fi
#    done
#    if [[ ${#fns[@]} -eq 0 ]]; then
#        printf "%s\n" "Dump file(s) do not exist. Exiting ..."
#        exit 1
#    fi
#    tag=()
#    if [[ $ns -eq 0 ]]; then
#        for (( i=0; i<${#fns[@]}; ++i )); do tag+=(${fns[$i]:np}); done
#    else
#        for (( i=0; i<${#fns[@]}; ++i )); do tag+=(${fns[$i]:np: -ns}); done
#    fi
#
#    IFS=$'\n'; tag=($(sort -n <<<"${tag[*]}")); unset IFS
#    fns=()
#    for each in ${tag[@]}; do fns+=("${prefix}"${each}"${suffix}"); done
#fi

#printf "%s\n" ${#fns[@]} > "${fn_fns}"
#for (( i=0; i<${#fns[@]}; ++i )); do
#    printf "%s\n" "${fns[$i]}" >> "${fn_fns}"
#done

#ibeg=1000; iend=1002; n=$(( iend-ibeg ))
##ibeg=1000; iend=${#fns[@]}; n=$(( iend-ibeg ))
#printf "%s\n" ${n} > "${fn_fns}"
#for (( i=$ibeg; i<$iend; ++i )); do
#    printf "%s\n" "${fns[$i]}" >> "${fn_fns}"
#done

#fn_fns=${dir}/traj/markers

./ligan \
    --job "zdf"  \
    --fn_cfg "${fn_cfg}"         \
    --fn_molinfo "${fn_molinfo}" \
    --fn_fns "${fn_fns}"         \
    --col_nums 1,2,3,4,5         \
    --coord_desc 'sw','sw','sw'  \
    --frame 1 --to -1 --stp 1  \
    --read_markers


#rm -f ${fn_fns}   

#{zdf, lan, lbm}
