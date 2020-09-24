#global parameters 
cp2k_exe=`which cp2k`
Tbegin=250
Tend=2500
Rate=20             # 1 epoch, temperature + 20 K
Step_per_epoch=20   # 1 epoch, 20 ionic steps

function cp2k-execute(){
    ${cp2k_exe} -i $1 -o cp2k.log
}

function cp2k_params() 
{
    #define parameters 
    step=$1
    timestep=1.0 
    temperature=$2
    pressure=$3 #atm
    freq_coord=$1
    freq_velocity=$1
    freq_restart=$1
    restart=$4

    # prepare input
    #atom->bar
    pressure=$(echo "scale=3; $pressure*1.01325"| bc)
    cp -f cp2k_initial.inp cp2k.inp;
    for i in step timestep temperature pressure freq_coord freq_velocity freq_restart;
        do 
            eval tmp=\$$i;
            sed -i "s/&{$i}/$tmp/g" cp2k.inp;
        done
        
    if [ ! $restart ];then
        true
    else
        echo  "&EXT_RESTART" >> cp2k.inp;
        echo  "RESTART_FILE_NAME $restart" >> cp2k.inp;
        echo  "&END EXT_RESTART" >> cp2k.inp;
    fi
}

#first run to obtain restart and wavefunction
mkdir -p $Tbegin;
cp2k_params 1 $Tbegin 0.0;
cd $Tbegin;
cp ../cp2k.inp ../coord.inc ./;
cp2k-execute cp2k.inp;
cp -f cp2k-1.restart ../cp2k_prev.restart;
cp cp2k-RESTART.wfn ../pre-cp2k-RESTART.wfn;
cd ../;
rm -rf $Tbegin;

#Doing epoches
total_epoches=$(echo "scale=0; ($Tend-$Tbegin)/$Rate"| bc)

for i in `seq 0 $total_epoches`;
    do 
        temp_now=$(echo "scale=3; $Tbegin+$Rate*$i"| bc);
        mkdir -p $temp_now;
        cp2k_params $Step_per_epoch $temp_now 0.0 cp2k_prev.restart;
        cd $temp_now;
        cp ../cp2k.inp ../coord.inc ../cp2k_prev.restart ./;
        mv ../pre-cp2k-RESTART.wfn ./cp2k-RESTART.wfn;
        cp2k-execute cp2k.inp;
        cp -f cp2k-1.restart ../cp2k_prev.restart;
        python ../cp2k2gro.py;
        cp system.gro ../${temp_now}.gro;
        cp cp2k-RESTART.wfn ../pre-cp2k-RESTART.wfn;
        cd ../;
        rm -rf $temp_now;
    done
 
rm -rf pre-cp2k-RESTART.wfn cp2k_prev.restart;

