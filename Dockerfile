FROM scevan_pyomics

RUN mkdir -p /home/f/feiler/dbenchSCEVAN
WORKDIR /home/f/feiler/dbenchSCEVAN
COPY . .

RUN yum install -y pip
RUN yum install -y python3-devel
RUN yum install -y libpng-devel
RUN pip install --no-cache-dir -r requirements.txt
RUN R -e "install.packages('Seurat',dependencies=TRUE, repos='http://cran.rstudio.com/')"

CMD ["python3", "/home/f/feiler/dbenchSCEVAN/run_scevan.py"]