FROM gcc:8.3 as gcc

RUN wget https://www.python.org/ftp/python/3.7.3/Python-3.7.3.tgz \
	&& tar -xvzf Python-3.7.3.tgz \
	&& cd Python-3.7.3 && ./configure && make && make install && cd ../ && rm -rf Python-3.7.3 && rm Python-3.7.3.tgz
RUN cd /usr/bin && unlink python && ln -s /usr/local/bin/python3 python && cd /
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_23.1.0-1-Linux-x86_64.sh \
	&& bash Miniconda3-py37_23.1.0-1-Linux-x86_64.sh -b -p /root/miniconda3 && /root/miniconda3/bin/conda init && cd / && rm Miniconda3-py37_23.1.0-1-Linux-x86_64.sh
ENV PATH /root/miniconda3/bin:$PATH

###set up environment
RUN apt update && apt-get install -y default-jre
RUN wget https://cran.r-project.org/src/base/R-3/R-3.6.1.tar.gz \
	&& tar -xvzf R-3.6.1.tar.gz && cd R-3.6.1 && ./configure --prefix=/usr/local && make && make install && cd / && rm -rf R-3.6.1.tar.gz R-3.6.1
RUN Rscript -e "install.packages(c('BiocManager'), repos='http://cran.us.r-project.org');BiocManager::install(c('CAGEr', 'BSgenome.Hsapiens.UCSC.hg38') );"

#####install ANGEL
RUN conda install -c anaconda gxx_linux-64 gcc_linux-64 cython
RUN conda install biopython=1.77 scikit-learn setuptools
RUN git clone https://github.com/Magdoll/ANGEL.git \
       && cd ANGEL \
       && python setup.py build \
       && python setup.py install
#install ANGEL end

RUN curl -LJO https://github.com/fulcrumgenomics/fgbio/releases/download/1.0.0/fgbio-1.0.0.jar

RUN conda install -c bioconda samtools bedtools
###RUN apt install less

##install kentUCSC 
RUN apt install -y rsync && rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/ /usr/local/bin 
##



RUN conda install -c bioconda htslib
RUN conda install -c bioconda seqtk
RUN conda install -c bioconda pysam
RUN conda install -c bioconda cd-hit

ENV PATH $PATH:/script






