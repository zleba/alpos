tag=AExt*

pwd=$PWD
cd $pwd/variants

for f in ${tag}.str
do
    $pwd/run.py $f
    sleep 0
done

#exit
#Submit them all

cd $pwd

for f in variants/${tag}_dir
do
   echo $f
   condor_submit  -batch-name `basename $f` dirName=$f  alpos.submit
done
