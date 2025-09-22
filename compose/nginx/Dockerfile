FROM nginx:1.19-alpine

RUN rm /etc/nginx/conf.d/default.conf
COPY ./compose/nginx/nginx.conf /etc/nginx/conf.d
COPY ./compose/nginx/maintenance.html /usr/share/nginx/html/maintenance.html