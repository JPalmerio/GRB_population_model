# ___________________________________________________________________________
#
# Jesse Palmerio - March 2020
#
# Makefile used to simplify work with docker images
# The base image has the necessary python packages and fortran libraries
# installed (to avoid rebuilding which takes a while)
# The core image has the LGRB population code installed on top of the
# base image 
# ___________________________________________________________________________

VERSION=0.1
REGISTRY=svomtest.svom.eu:5543
BASE_IMG=LGRBpop_base:${VERSION}
CORE_IMG=LGRBpop_core:${VERSION}
FULL_IMG=LGRBpop:${VERSION}
CONT_NAME=LGRBpop_container
STACK_NAME=LGRBpop

.PHONY: base_img core bashrun deploy clean log


base_img:
    @echo '--'
    @echo 'Docker build base image from base_image_Dockerfile'
    @echo '--'
    @docker build --tag=${BASE_IMG} -f base_img/base_img_Dockerfile .
    @afplay /System/Library/Sounds/Morse.aiff -v 10

core:
    @echo '--'
    @echo 'Docker build core image without cache from Dockerfile ${PWD}'
    @echo '--'
    @docker build --no-cache --rm --tag=${CORE_IMG} --build-arg BASE_IMAGE=${BASE_IMG} .

bashrun:
    @echo '--'
    @echo 'Docker run ${CORE_IMG} in container: ${CONT_NAME}'
    @echo '--'
    @docker run --rm -it --entrypoint=ash --name=${CONT_NAME} ${CORE_IMG}


deploy:
    @echo '--'
    @echo 'Docker stack deployment'
    @echo '--'
    @docker stack deploy --with-registry-auth -c configuration.yml --prune ${STACK_NAME}

clean:
    @if [ $$(docker service ls -qf name=${STACK_NAME}|wc -l) -gt 0 ]; then\
        echo '--';\
        echo 'Old stack cleansing';\
        echo '--';\
        docker stack rm ${STACK_NAME};\
    fi
    @while [ $$(docker ps -qf name=${STACK_NAME}|wc -l) -gt 0 ]; do\
        sleep 1;\
    done
    @sleep 2

log:
    @echo '--'
    @echo 'Service log'
    @echo '--'
    @docker service logs ${STACK_NAME}_pipeline -f --raw


notrunc:
    @echo '--'
    @echo 'No trunc'
    @echo '--'
    @docker service ps ${STACK_NAME}_pipeline --no-trunc

coverage:
    @echo '--'
    @echo 'Running coverage'
    @echo '--'
    coverage run --source src/ -m pytest
    @coverage xml -o reports/coverage.xml


#______________________________________________________________________________
