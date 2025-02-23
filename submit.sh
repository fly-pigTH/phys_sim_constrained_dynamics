# Create the submission folder.
mkdir -p submission

# Copy the source code.
cp video.mp4 submission/
RET=$?
if [ $RET != 0 ]; then
    echo "$(tput setaf 1)Submission failed. Please check the output."
    exit
fi

cp cpp/link/src/link.cpp submission/
cp cpp/joint/src/binary_hinge_joint.cpp submission/
cp cpp/sim/src/simulator_step.cpp submission/
cp cpp/contact/src/collision_detector.cpp submission/

zip -r submission.zip submission