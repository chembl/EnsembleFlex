
# Container for building the environment

# Installation of pip dependencies from (explicit) lockfiles is only supported by the 'conda-lock install' and
# 'micromamba install' commands (not conda or mamba). Mamba is faster than conda.
# So, we use a micromamba base image. Documentation on the micromamba docker images can be found here:
# https://micromamba-docker.readthedocs.io/en/latest/quick_start.html
FROM mambaorg/micromamba:1.5.8

ENV DEBIAN_FRONTEND=noninteractive

# This fix: libGL error: No matching fbConfigs or visuals found
ENV LIBGL_ALWAYS_INDIRECT=1

# set the working directory
WORKDIR /app

# install conda environment from lock file
COPY --chown=$MAMBA_USER:$MAMBA_USER conda-lock.yml /tmp/conda-lock.yml
# installs into the base environment
RUN micromamba install --name base --yes --file /tmp/conda-lock.yml \
    && micromamba clean --all --yes

ENV PATH="/env/bin:${PATH}"

# copy scripts to the working directory
# (as last steps so that when there are changes to the scripts the above steps don't need to be repeated as they are cached)
COPY . /app
# copy configuration files from .streamlit folder (e.g. config.toml)
COPY ./src/streamlit_app/.streamlit/ /app/.streamlit/

# expose Streamlit’s (default) network port
EXPOSE 8501

# tell Docker how to test a container to check that it is still working.
# (Your container needs to listen to Streamlit’s (default) port 8501.)
HEALTHCHECK CMD curl --fail http://localhost:8501/_stcore/health

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "streamlit", "run", "src/streamlit_app/streamlit_app.py", "--server.port=8501", "--server.address=0.0.0.0"]
