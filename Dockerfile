FROM ubuntu:20.04
RUN apt-get update
RUN apt-get install -y apt-utils
RUN apt-get install -y debconf
RUN apt-get install -y tzdata
RUN dpkg-reconfigure --frontend noninteractive tzdata
RUN apt-get install -y samtools bowtie2 bedtools
RUN apt-get install -y python3-pip
RUN pip3 install pandas biopython pysam
ADD . .
