#
# docker build -t registry.nersc.gov/m4616/raknop/nersc-desi-gaia-dr3-server .
#

FROM rknop/devuan-daedalus-rknop AS base

MAINTAINER Rob Knop <raknop@lbl.gov>

SHELL [ "/bin/bash", "-c" ]

RUN apt-get update \
    && DEBIAN_FRONTEND="noninteractive" apt-get -y upgrade \
    && DEBIAN_FRONTEND="noninteractive" TZ="US/Pacific" apt-get -y install -y \
         python3 python3-venv \
    && apt-get -y autoremove \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ======================================================================
# pip installs a full dev environment, which we don't want
#  in our final image.  (400 unnecessary MB.)

FROM base AS build

RUN apt-get update \
    && DEBIAN_FRONTEND="noninteractive" apt-get install -y python3-pip

RUN mkdir /venv
RUN python3 -mvenv /venv

RUN source /venv/bin/activate \
  && pip install \
       gunicorn flask numpy astropy healpy

RUN mkdir /tmp/build
RUN mkdir /code

COPY . /tmp/build
WORKDIR /tmp/build

RUN make INSTALLDIR=/code install

# ======================================================================

FROM base AS final

COPY --from=build /venv/ /venv/
ENV PATH=/venv/bin:$PATH

COPY --from=build /code/ /code/
WORKDIR /code

# This next one gets bind mounted to /global/cfs/cdirs/desi/target/gaia_dr3/healpix
RUN mkdir /gaia

CMD [ "gunicorn", "-w", "4", "-b", "0.0.0.0:8080", "--timeout", "0", \
      "webservice:app" ]
