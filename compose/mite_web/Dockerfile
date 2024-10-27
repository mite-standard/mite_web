FROM python:3.12.0b1-slim-buster

ENV PYTHONUNBUFFERED 1

ENV INSTALL_PATH /mite_web
RUN mkdir -p $INSTALL_PATH
WORKDIR $INSTALL_PATH

COPY ./mite_web .

RUN pip install --upgrade pip
RUN pip install .

RUN chmod +x ./entrypoint_docker.sh
RUN chmod +x ./prepare_mite_data.py
RUN python3 ./prepare_mite_data.py

EXPOSE 8004

CMD ["./entrypoint_docker.sh"]
