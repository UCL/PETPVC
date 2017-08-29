# PETPVC

FROM alpine:3.6

MAINTAINER Benjamin A. Thomas <b.a.thomas@ucl.ac.uk>
LABEL Description="PETPVC" 

RUN apk --update add ca-certificates && apk upgrade

RUN apk add alpine-sdk cmake && \
    mkdir -p /opt/ITK/BUILD && cd /opt/ITK && \
    git clone -b 'release' --single-branch --depth=1 https://github.com/InsightSoftwareConsortium/ITK.git && \
    cd /opt/ITK/BUILD && \
    cmake -DBUILD_TESTING:BOOL=OFF -DModule_ITKReview:BOOL=ON -DCMAKE_BUILD_TYPE:STRING=Release /opt/ITK/ITK && \
    make && \
    mkdir -p /opt/PETPVC/BUILD && cd /opt/PETPVC && \
    git clone https://github.com/UCL/PETPVC && \
    cd BUILD && cmake -DCMAKE_BUILD_TYPE:STRING=Release -DITK_DIR=/opt/ITK/BUILD /opt/PETPVC/PETPVC && \
    make && make install && ctest && cd ~ && \
    rm -rf /opt/PETPVC && rm -rf /opt/ITK && \
    apk del --purge alpine-sdk cmake && apk add libstdc++
    
