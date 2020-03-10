ARG BASE_IMAGE=null
FROM ${BASE_IMAGE}

MAINTAINER Jesse PALMERIO <jesse.palmerio@obspm.fr>

ENV PATH /home/src/:${PATH}
ENV PYTHONPATH /home/src/:${PYTHONPATH}


COPY ./src/ /home/src
COPY ./data/ /home/data
COPY ./init/ /home/init
COPY ./observational_constraints/ /home/observational_constraints

RUN mkdir /home/model_outputs
VOLUME /home/model_outputs

WORKDIR /home/
RUN cd src; make; cd ..

ENTRYPOINT ["python3","src/simple_example.py"]