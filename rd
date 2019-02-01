#/bin/sh
image=zleba/rootlhapdf:system

rel=$(python -c "import os.path; print os.path.relpath('$PWD', '$PROJECT_DIR')")

if [ ! -z "$DISPLAY" ]
then
    #sudo docker run -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix -v $PWD:/excercise --rm -it --user $(id -u) $image /bin/bash -c "cd /excercise/; $*"
    
else
    sudo docker run  -v  $PROJECT_DIR:/excercise --rm -it --user $(id -u) $image /bin/bash -c "cd /excercise/;source setup.sh;cd $rel; $*"
fi
