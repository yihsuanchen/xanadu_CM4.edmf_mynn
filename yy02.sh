#!/bin/bash 
#===================================
#  Bourne-Again shell script
#
#  Description:
#
#
#  Usage:
#
#
#  History:
#
#  Author:
#    Yi-Hsuan Chen
#    yihsuan@umich.edu
#===================================

#************************
#  Description:
#    create commands to execute NCL files with reading inputs from command lines
#************************

###################
# user setting
###################

### ncl files ###
nclfiles="y2-plot-profile-AMIP_i27_j01_IndOcn.ncl"
#nclfiles=(\
#          ".ncl" \
#          ".ncl" \
#         )

### input variable names ###
#vars_input=("infilename_noScat" "infilename_Scat" "infilename_ScatFIR" "plotname_step" "plotname_file" "varvars")
#vars_input=("plotname_step" "plotname_file" "varvars" "filehead_Scat" "filehead_noScat")
vars_input=("option_var")

### group input variables, only support g01~g09 ###
vars_group=("g01")
#vars_group=("g01" "g02" "g04" "g03" "g03")
#vars_input=("thl" "qall" "tke" "tdt" "qdt" "udt" "vdt")

#*** set choice of each input variables ***
function read_choices_input {
  local var_in=$1

  if [ $var_in == "lskjdflksjdlfkjsaldflajskdfl" ]; then
    choice_out=""
  elif [ $var_in == "option_var" ]; then
    choice_out=("thl" "qall" "tke" "tdt" "qdt" "udt" "vdt")
  elif [ $var_in == "plotname_step" ]; then
    choice_out=("ANN" "DJF" "MAM" "JJA" "SON")
  elif [ $var_in == "plotname_file" ]; then
    choice_out=("global")
    #choice_out=("global" "Antarctic" "midlats_SH" "tropics" "midlats_NH" "Arctic")
    #choice_out='(/q1,q2/)'
  elif [ $var_in == "varvars" ]; then
#    choice_out=(\
#               '(/"FLDS","FLDSC","FLUS","FSDS","FSDSC","FSNS","FSNSC","FSUS","FLUT","FLUTC","FSNTOA","FSNTOAC"/)' \
#               '(/"SHFLX","LHFLX","TS","TREFHT","ICEFRAC","RADSURF"/)' \
#               '(/"CLDLOW","CLDMED","CLDHGH","CLDTOT","LWCF","SWCF","LWCF_SURF","SWCF_SURF"/)' \
#               '(/"PRECT","PRECC","PRECL","EVAP","P_minus_E","TMQ","LWCOOL_TOA_TO_SURF","SWWARM_TOA_TO_SURF","RADCOOL_TOA_TO_SURF"/)' \
#                )
    choice_out=(\
               '(/"TREFHT","FLDS"/)' \
                )
  elif [ $var_in == "filehead_Scat" ]; then
    #choice_out="../data/h01-cesm111-FC5-mc6_rtr2_Scat.cam.h0.yy06_35-climo_"  # will be set later
    choice_out="../data/c10-cesm111-E2000_rrtmg_mc6_rtr2.cam.h0.yy06_35-climo_"  # will be set later
  elif [ $var_in == "filehead_noScat" ]; then
    #choice_out="../data/h02-cesm111-FC5-mc6_rtr2_noScat.cam.h0.yy06_35-climo_"  # will be set later
    choice_out="../data/c11-cesm111-E2000_rrtmg_mc6_rtr2_noScat.cam.h0.yy06_35-climo_"  # will be set later
  else
    echo "ERROR: input variable [$var_in] does not exist"
    echo "program stop"
    exit 1
  fi

  echo "${choice_out[@]}!${#choice_out[@]}"
} # end read_choices_input

#**********************************************
#  Example
#
#    vars_input=("vv1" "vv2" "vv3")
#    vars_group=("g01" "g01" "g02")   # vv1 & vv2 in a group
#    vv1=("ANN"  "DJF"  "MAM"  "JJA"  "SON")
#    vv2=("ANN0" "DJF0" "MAM0" "JJA0" "SON0")
#    vv3=("1" "2" "3")
#
#    NCL commands
#      ncl vv1=ANN vv2=ANN0 vv3=1
#      ncl vv1=ANN vv2=ANN0 vv3=2
#      ncl vv1=ANN vv2=ANN0 vv3=3
#      ncl vv1=DJF vv2=DJF0 vv3=1
#      ....
#**********************************************

##################
# program start
##################
#set -x

temp=`date +%Y%m%d%H%M%S`

nvar_input=${#vars_input[@]}
nvar_group=${#vars_group[@]}
nnclfiles=${#nclfiles[@]}
ftemp01="./sskkccoo.$temp.01"

cat >> $ftemp01 << EOF1
#!/bin/bash 
set -x

EOF1

#---------------------------
# check nclfiles with user
#---------------------------

echo "------------------------"
echo "Execute these NCL files"
echo ""
for ((i=0; i<$nnclfiles; i=i+1))
do
  j=$(($i+1))
  file1=${nclfiles[$i]}
  if [ -f $file1 ]; then 
    echo "  input files        : $j/$num_infile, [${nclfiles[$i]}]"
  else
    echo "ERROR: file [$file1] does not exist"
    echo "Program stop"
    exit 2
  fi
done
echo "------------------------"

#--------------------------
# read variables & groups
#--------------------------

#set -x

if [ $nvar_input -ne $nvar_group ]; then
  echo "ERROR: #vars_input is not equal to #vars_group"
  echo "  vars_input: #${nvar_input}, [${vars_input[@]}]"
  echo "  vars_input: #${nvar_group}, ${vars_group[@]}"
  exit 1
fi

# initialize number of group
g01=1; g02=1; g03=1; g04=1; g05=1
g06=1; g07=1; g08=1; g09=1

for ((vv=0; vv<$nvar_input; vv=vv+1))
do
  var1=${vars_input[$vv]}
  group1=${vars_group[$vv]}

  ppp=$(read_choices_input $var1)
  choice1=`echo $ppp | cut -d '!' -f 1` 
  num1=`echo $ppp | cut -d '!' -f 2` 
  #echo $var1, $num1, $ppp

  # check number of group is consistent
  if [ $group1 == "g01" -a $g01 -ne 1 -a $g01 -ne $num1 ] ; then exit_code="201" ; v00=$v01 ; g00=$g01 ; fi
  if [ $group1 == "g02" -a $g02 -ne 1 -a $g02 -ne $num1 ] ; then exit_code="201" ; v00=$v02 ; g00=$g02 ; fi
  if [ $group1 == "g03" -a $g03 -ne 1 -a $g03 -ne $num1 ] ; then exit_code="201" ; v00=$v03 ; g00=$g03 ; fi
  if [ $group1 == "g04" -a $g04 -ne 1 -a $g04 -ne $num1 ] ; then exit_code="201" ; v00=$v04 ; g00=$g04 ; fi
  if [ $group1 == "g05" -a $g05 -ne 1 -a $g05 -ne $num1 ] ; then exit_code="201" ; v00=$v05 ; g00=$g05 ; fi
  if [ $group1 == "g06" -a $g06 -ne 1 -a $g06 -ne $num1 ] ; then exit_code="201" ; v00=$v06 ; g00=$g06 ; fi
  if [ $group1 == "g07" -a $g07 -ne 1 -a $g07 -ne $num1 ] ; then exit_code="201" ; v00=$v07 ; g00=$g07 ; fi
  if [ $group1 == "g08" -a $g08 -ne 1 -a $g08 -ne $num1 ] ; then exit_code="201" ; v00=$v08 ; g00=$g08 ; fi
  if [ $group1 == "g09" -a $g09 -ne 1 -a $g09 -ne $num1 ] ; then exit_code="201" ; v00=$v08 ; g00=$g09 ; fi

  if [ $exit_code -a $exit_code == "201" ]; then
    echo "ERROR: group [$group1] has different # of input choices"
    echo "  var=[$v00, $var1], #choice=[$g00, $num1]"
    echo "check function [read_choices_input] to modify"
    echo ""
    exit 1
  fi

  if [ $group1 == "g01" ] ; then g01=$num1 ; v01=$var1 ; fi
  if [ $group1 == "g02" ] ; then g02=$num1 ; v02=$var1 ; fi
  if [ $group1 == "g03" ] ; then g03=$num1 ; v03=$var1 ; fi
  if [ $group1 == "g04" ] ; then g04=$num1 ; v04=$var1 ; fi
  if [ $group1 == "g05" ] ; then g05=$num1 ; v05=$var1 ; fi
  if [ $group1 == "g06" ] ; then g06=$num1 ; v06=$var1 ; fi
  if [ $group1 == "g07" ] ; then g07=$num1 ; v07=$var1 ; fi
  if [ $group1 == "g08" ] ; then g08=$num1 ; v08=$var1 ; fi
  if [ $group1 == "g09" ] ; then g09=$num1 ; v09=$var1 ; fi

  vars_choice[$vv]=$choice1
  num_choice[$vv]=$num1
done

# execute each nclfile
for file1 in ${nclfiles[@]}
do

  for ((ii01=1; ii01<=$g01; ii01=ii01+1)) ; do
  for ((ii02=1; ii02<=$g02; ii02=ii02+1)) ; do 
  for ((ii03=1; ii03<=$g03; ii03=ii03+1)) ; do 
  for ((ii04=1; ii04<=$g04; ii04=ii04+1)) ; do 
  for ((ii05=1; ii05<=$g05; ii05=ii05+1)) ; do 
  for ((ii06=1; ii06<=$g06; ii06=ii06+1)) ; do 
  for ((ii07=1; ii07<=$g07; ii07=ii07+1)) ; do 
  for ((ii08=1; ii08<=$g08; ii08=ii08+1)) ; do 
  for ((ii09=1; ii09<=$g09; ii09=ii09+1)) ; do 

     print_command=""
     exe_command=""
     for ((vv=0; vv<$nvar_input; vv=vv+1))
     do
       var1=${vars_input[$vv]}
       group1=${vars_group[$vv]}
       choice_all=${vars_choice[$vv]}

       if [ $group1 == "g01" ] ; then gg0=$ii01 ; fi
       if [ $group1 == "g02" ] ; then gg0=$ii02 ; fi
       if [ $group1 == "g03" ] ; then gg0=$ii03 ; fi
       if [ $group1 == "g04" ] ; then gg0=$ii04 ; fi
       if [ $group1 == "g05" ] ; then gg0=$ii05 ; fi
       if [ $group1 == "g06" ] ; then gg0=$ii06 ; fi
       if [ $group1 == "g07" ] ; then gg0=$ii07 ; fi
       if [ $group1 == "g08" ] ; then gg0=$ii08 ; fi
       if [ $group1 == "g09" ] ; then gg0=$ii09 ; fi

       choice1=`echo $choice_all | cut -d ' ' -f $gg0`
       print_command="${print_command} '$var1=\"$choice1\"'"
       exe_command="${exe_command} $var1=\"$choice1\""

       print_command=`echo $print_command | sed s,\"\(,\(,g | sed s,\)\",\),g`
       exe_command=`echo $exe_command | sed s,\"\(,\(,g | sed s,\)\",\),g`
     done  # end loop of nvar_input

     # save ncl commands
     echo "ncl $print_command $file1 || exit 5" >> $ftemp01

     #echo $print_command
     #echo $exe_command
  
     #if [ $opt_run -a $opt_run == "print" ]; then
     #  echo "ncl $print_command $file1"
     #else
       #echo ""
       #echo "ncl $print_command $file1" >>
       #echo "ncl $exe_command $file1" >> $ftemp01
       #echo ""
       #ncl $exe_command $file1 || exit 5
     #fi

  done ; done ; done ; done ; done
  done ; done ; done ; done

done  # end loop of nclfiles

#-----------------------
# execute NCL commands
#-----------------------
  cat $ftemp01

  echo " "
  read -p "Is is correct? (y/n)  " choice
  echo " "
  if [ $choice == "y" ]; then
    chmod 755 $ftemp01 || exit 1
    $ftemp01 || exit 2
  else
    echo "Cancel by user"
    echo "program stop"
  fi

  rm $ftemp01 || exit 3

exit 0

