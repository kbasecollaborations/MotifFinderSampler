#!/bin/bash

dest=/kb/deployment/bin/meme
curl -O http://meme-suite.org/meme-software/5.0.1/meme_5.0.1_1.tar.gz
tar zxf meme_5.0.1_1.tar.gz
cd meme-5.0.1
./configure --prefix=$dest
make
make install
export PATH=$dest/bin:$PATH
cd ..
