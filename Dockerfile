#
# docker build -t registry.nersc.gov/m4616/raknop/nersc-desi-gaia-dr3-server:yyyymmdd .
#

FROM debian:trixie-20260112 AS base
LABEL maintainer="Rob Knop <raknop@lbl.gov>"

SHELL [ "/bin/bash", "-c" ]

RUN apt-get update \
    && DEBIAN_FRONTEND="noninteractive" apt-get -y upgrade \
    && DEBIAN_FRONTEND="noninteractive" TZ="US/Pacific" apt-get -y install -y \
         net-tools procps python3 python3-venv \
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
       gunicorn flask numpy astropy healpy gevent

RUN mkdir /tmp/build
RUN mkdir /code

COPY . /tmp/build
WORKDIR /tmp/build

RUN make INSTALLDIR=/code install

# ======================================================================

FROM base AS final

COPY --from=build /venv/ /venv/
ENV PATH=/venv/bin:$PATH
# ...gunicorn seems to want to write a dotfile to my home directory?????????
# Ah, it's a control socket.  Why it writes it to $HOME, I don't know.
ENV HOME=/tmp

COPY --from=build /code/ /code/
WORKDIR /code

# This next one gets bind mounted to /global/cfs/cdirs/cosmo/data/gaia/dr3/healpix
RUN mkdir /data

CMD [ "gunicorn", "-w", "4", "-b", "0.0.0.0:8080", "-k", "gevent", "--timeout", "30", "webservice:app" ]
