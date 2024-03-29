#!/bin/bash

# Adapted from https://github.com/MaayanLab/Zika-RNAseq-Pipeline/blob/master/notebook.sh
# Strict mode
set -euo pipefail
IFS=$'\n\t'

# Set configuration defaults
: ${PASSWORD:=""}
: ${PEM_FILE:="/key.pem"}
: ${USE_HTTP:=0}

# Create the hash to pass to the IPython notebook, but don't export it so it doesn't appear
# as an environment variable within IPython kernels themselves
HASH=$(python -c "from IPython.lib import passwd; print(passwd('${PASSWORD}'))")
unset PASSWORD

if [ $USE_HTTP -ne 0 ]; then
  CERTFILE_OPTION=""
else
  # Create a self signed certificate for the user if one doesn't exist
  if [ ! -f $PEM_FILE ]; then
    openssl req -x509 -nodes -days 365 -newkey rsa:2048 -keyout $PEM_FILE -out $PEM_FILE \
      -subj "/C=XX/ST=XX/L=XX/O=dockergenerated/CN=dockergenerated"
  fi

  CERTFILE_OPTION="--certfile=$PEM_FILE"
fi

jupyter notebook \
	--port 8888 \
	--ip=0.0.0.0 \
        --NotebookApp.base_url='mirpipe' \
	$CERTFILE_OPTION \
	--NotebookApp.password="$HASH" \
        --allow-root

unset HASH
