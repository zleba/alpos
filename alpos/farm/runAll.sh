tag=GExt_nnlo_heraCfps4D*

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
   #echo $f
   nFiles=`ls  $f/steering.str?  $f/steering.str??  | wc -l`
   condor_submit  -batch-name `basename $f` dirName=$f nFiles=$nFiles  alpos.submit
done
