#!/bin/bash

set -xe

name=$(date +%s)
mkdir -p .backups
mv -i store.sqlite .backups/$name.sqlite
