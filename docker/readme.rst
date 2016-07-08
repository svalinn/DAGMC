How to Update the CI Docker Image
------------
Make changes to the Dockerfile as needed, then rebuild (assuming you are in thedirectory wherer the Dockfile is

    sudo docker build . 

Docker will now build the image according to the instructions in the file. This must now be updated, check the last built image by running

    sudo docker images

The last one is the one we want to push, so tag it

    sudo docker tag <id num> makeclean/dagmc-ci:latest

Set your Docker Hub username and email address

    sudo docker login --username=makeclean --email=andrew.davis@wisc.edu

Now push the file to DockerHub

    sudo docker push makeclean/dagmc-ci

CI is now updated, and should pull from here
